/** 
 * @package stringanimation
 */

define(['d3', 'knockout'], function (d3, ko) {
    
  return {
    dataset: ko.observable(),
    width: 0,
    height: 110,
    letterWidth: 1,

    /**
     * Duration of a single step of the animation in ms. Most animations take up 3 steps.
     * @property transitionDuration
     */ 
    transitionDuration: 200,
    
    /**
     * Duration of a scroll in ms
     * @property scrollDuration
     */ 
    scrollDuration: 500,
    
    /**
     * Determine the width of a single character by creating a hidden canvas.
     * @method charWidth
     * @param {String} font The font used in the animated string
     * @return {Integer} The width of a single letter in px
     */
    charWidth: function(font) {      
      var canvas = document.getElementById('graph-string-measurement-canvas');
      var context = canvas.getContext("2d");
      context.font = font;
      var metrics = context.measureText('m');
      return metrics.width;
    },
    
		/**
     * Transform the nested mutations with pieces to a flat block structure in which 
     * also the unchanged parts of the DNA-string are included.
     * @method transformData
     * @return {Array} Transformed data
     */
    transformData: function () {
      var self = this;
      
      var ref = this.dataset().reference_sequence();

      // Cut sample string in blocks
      var sample_sequence = this.dataset().sample_sequence();
      var reference_sequence = this.dataset().reference_sequence();
      var last_end_position = 1;
      self.blocks = [];
      self.block_index = {};
      var id_index = 0;
      
      // Iterate over all mutations and add them to the array.
      self.dataset().allele_description().forEach(function(mutation) {
        // make the string between last known end position and start of new mutation if > 0
        var unmodified_length = mutation.reference_start() - last_end_position;
        
        // If there is a piece of unmodified text between the new and the last mutation, add this to the array
        if (unmodified_length > 0) {
          self.blocks.push({
            id: id_index,
            type: 'text',
            mutation_type: null,
            sequence_before: reference_sequence.substr(last_end_position-1, unmodified_length),
            sequence_after: reference_sequence.substr(last_end_position-1, unmodified_length),
            state: 'before',
            position: 0,
            position_origin: null,
            _mutation: undefined,
            _piece: undefined,
          });
          id_index++;
        }
        
        var mutation_sample_length = mutation.sample_end() - mutation.sample_start();
        
        // Add all pieces to the block array
        mutation.pieces()
          .forEach(function(piece) {
            self.blocks.push({
              id: id_index,
              type: piece.type(),
              mutation_type: mutation.type(),
              sequence_before: (piece.type() === 'del') ? piece.sequence() : '',
              sequence_after: (piece.type() === 'del') ? '' : piece.sequence(),
              state: 'before',
              position: 0,
              position_origin: piece.origin_start(),
              _mutation: mutation, // maintain reference to the mutation
              _piece: piece, // maintain reference to the piece
            });
            id_index++;
            piece._block = self.blocks[self.blocks.length-1];
          });
                  
        // Due to differences in position numbering for different mutation types, this statement is needed
        if (mutation.type() == 'ins' || mutation.type() == 'dup') {
          last_end_position = mutation.reference_end();
        } else {
          last_end_position = mutation.reference_end();
        }
      }); // end forEach mutation
      
      // After the last mutation, there still is one last piece of text.
      if (last_end_position < reference_sequence.length) {
        self.blocks.push({
          id: id_index,
          type: 'text',
          mutation_type: null,
          sequence_before: reference_sequence.substr(last_end_position-1),
          sequence_after: reference_sequence.substr(last_end_position-1),
          state: 'before',
          position: 0,
          position_origin: null,
            _mutation: undefined,
            _piece: undefined,          
        });
        id_index++;
      }
    },
    
    /**
     * Called by Durandal upon activating the ViewController
     * @method activate
     */
		activate: function() {
			var self = this;
			self.dataset(ko.utils.unwrapObservable(require('viewmodels/visualisation').dataset));
		},
		
    /**
     * Called by Durandal upon attachment of this view to the DOM
     * @method attached
     */
		attached: function() {
			var self = this;
      self.letterWidth = self.charWidth('30px Courier New');
			
			self.render();
		},
		
    /**
     * Reset the string to the original start position
     * @method resetString
     * @param {string} font The font used in the animated string
     */
		resetString: function() {
      var self = this;
      self.blocks.forEach(function (block) { block.state = 'before'; });
      self.render();
      self.recalculatePositions();
      self.rerender();
		},
		
    /**
     * Draw the non-changing parts of the animation: both the reference and sample string and add all elements
     * @method render
     */
    render: function() {
			var self = this;
			var ref = this.dataset().reference_sequence();

      $('.graph-string-animation svg').empty();
			self.svg = d3.select('.graph-string-animation svg');

      // set size of svg element to size of longest string length
      self.svg.attr('width', Math.max(self.dataset().reference_sequence().length, self.dataset().sample_sequence().length) * self.letterWidth);
      
			self.transformData();
      
      // Draw all the blocks
			var blocks = self.svg.selectAll('text.block').data(self.blocks, function(d) { return d.id; });
      blocks.enter().append('text').attr('class',function(d) { return 'state-before block ' + d.type + ' mutation-' + d.mutation_type;});
      blocks
        .attr('x', function(d) { return d.position * self.letterWidth; })
        .attr('y', 65)
				.text(function(d) { return d.sequence_before; });

			// temporary ruler
			self.svg.selectAll('text.ruler')
        .data(this.dataset().reference_sequence().split(''), function(d) { return d.id; })
        .enter().append('text').attr('y', 16).attr('x', function(d, i) { return self.letterWidth * (i-1)+6; }).text(function(d,i) {  return ((i) % 5 === 0) ? i : ''; })
        .attr('class', 'ruler');
		
      self.svg.append('text').attr('x', 4).attr('y', 36).attr('class', 'reference-string').text(this.dataset().reference_sequence().split('').join(' '));
      self.svg.append('text').attr('x', 4).attr('y', 85).attr('class', 'sample-string').text(this.dataset().sample_sequence().split('').join(' '));
    },

    /**
     * Calculate the positions of the blocks using string lengths and state of each block.
     * @method recalculatePositions
     */
    recalculatePositions: function () {
      var self = this;
      var position = 0;
      self.blocks.forEach( function(block) {
        block.position = position;
        
        // Special case to fix 'closing' the gap in animation of substitutions and invertions
        if (block.mutation_type === 'subst' || block.mutation_type === 'inv') {
          position += block.sequence_after.length;
        } else {
          position += (block.state == 'before') ? block.sequence_before.length : block.sequence_after.length;
        }
      });
    },

    /** 
     * Start animation of a single piece. Returns a promise.
     * - find piece in block array
     * - start scrolling the view to the correct position
     * - change state of element to 'after'
     * - recalculate positions
     * - refresh data display
     * 
     * @method animate
     * @param {Object} piece The piece (part of mutation) that should be animated
     * @return {Promise} A Promise that resolves when animation is complete
     * @async
     **/
    animate: function(piece) {
      var self = this;
      
      // Find the corresponding block
      var block = piece._block;
      
      return self.scrollToOffset(block.position)
        .then(function() {
          block.state = 'after';
          self.recalculatePositions();
          return self.rerender();
        });
    },
    
    /**
     * Scroll the view to a certain offset
     * @method scrollToOffset
     * @param {Integer} offset String position offset to scroll to
     * @return {Promise} A Promise that is resolved when scolling is complete
     * @async
     */
    scrollToOffset: function(offset) {
      var dfd = $.Deferred();
      $('.graph-string-animation').animate({
          scrollLeft: this.letterWidth * offset - 120
        }, {
          duration: this.scrollDuration, 
          easing: 'swing',
          done: function() {
            dfd.resolve();
          }
        });
      return dfd.promise();
    },
    
    /**
     * Rerender the string, only updating the updateble parts of the animation.
     * Animates blocks that have a state change
     * @method rerender
     * @return {Promise} A promise that resolves when all animations are complete
     * @async
     */
    rerender: function () {
      var self = this;
        
      // The y positions of string elements.
      var y = 65;
      var y_below = 100;
      var y_top = 30;
      
      // Function that calls a callback at the end of all animations within a d3 chain
      // @function endall
      var endall = function(transition, callback) { 
        var n = 0;
        transition
          .each(function() { ++n; }) 
          .each("end", function() { if (!--n) callback.apply(); }); 
      };
      
      // Select all blocks in the view
      var blocks = self.svg.selectAll('text.block').data(self.blocks, function(d,i) { return d.id; });
      
      // Select blocks that will be changed
      var transitions = self.svg.selectAll('text.block.state-before').data(self.blocks, function(d,i) { return d.id; })
        .filter(function(d, i) { return d.state === 'after'; });
      

      //  Select all deleted blocks that are changing state
      var dfd_del = $.Deferred();
      var del_transitions = transitions.filter(function(d) { return d.type === 'del'; });
      
      // If there are no delete transitions, just resolve the promise
      if (del_transitions.empty()) {
        dfd_del.resolve();
      } else {
        del_transitions
				  .transition().duration( self.transitionDuration )
					  .attr('y', y_below)
          .transition().duration( self.transitionDuration )
					  .style('fill-opacity', 1e-6)
            .text(function(d) { return d.sequence_after; })
            .attr('class', function(d) { return 'block state-after ' + d.type; })
          .transition().duration( self.transitionDuration )
            .call(endall, dfd_del.resolve);

        blocks.filter(function(d) { return d.type === 'text'; })
          .transition().duration( self.transitionDuration * 2 )
          .transition().duration( self.transitionDuration )
            .attr('x', function(d) { return d.position * self.letterWidth; });
      }
       
      // Transition for inserts
      var dfd_ins = $.Deferred();
      var ins_transitions = transitions.filter(function(d) { return d.type == 'ins'; });
      if (ins_transitions.empty()) {
        dfd_ins.resolve();
      } else {
        ins_transitions
            .attr('y', y_top)
            .attr('x', function(d) { return d.position * self.letterWidth; })
            .style('fill-opacity', 1e-6)
            .text(function(d) { return d.sequence_after; })
          .transition().duration(self.transitionDuration)
            .style('fill-opacity', 1)
          .transition().duration(self.transitionDuration)
            .attr('y', y)
          .transition().duration( self.transitionDuration)
            .call(endall, dfd_ins.resolve);

        blocks.filter(function(d,i) { return d.state === 'before'; })
          .transition().duration( self.transitionDuration )
          .attr('x', function(d) { return d.position * self.letterWidth; });
      }
                    
      // transition for transpositions
      var dfd_trans = $.Deferred();
      var trans_transitions = transitions.filter(function(d) { return d.type == 'trans'; });
      if (trans_transitions.empty()) {
        dfd_trans.resolve();
      } else {
        trans_transitions
            .attr('y', y_top)
            .attr('x', function(d) { return d.position_origin * self.letterWidth; })
            .style('fill-opacity', 1e-6)
            .text(function(d) { return d.sequence_after; })
          .transition().duration( self.transitionDuration)
            .style('fill-opacity', 1)
          .transition().duration( self.transitionDuration)
            .attr('x', function(d) { return d.position * self.letterWidth; })
          .transition().duration( self.transitionDuration)
            .attr('y', y)
          .transition().duration( self.transitionDuration)
          .call(endall, dfd_trans.resolve);

        blocks.filter(function(d,i) { return d.state === 'before'; })
          .transition().duration( self.transitionDuration )
          .attr('x', function(d) { return d.position * self.letterWidth; });
      }
      
      // After all animations are complete, update the state and x-position of entire string
      return $.when(dfd_del, dfd_trans, dfd_ins).then(function() {
        blocks.attr('x', function(d) { return d.position * self.letterWidth; });
        blocks.attr('class', function(d) { return 'block ' + 'state-'+d.state + ' ' + d.type + ' mutation-' + d.mutation_type; });
      }); 
    },

  };
});
