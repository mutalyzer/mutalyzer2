define(['d3', 'durandal/app', 'knockout'], function (d3, app, ko) {

  return {
    root: null,
    dataset: ko.observable(),
    width: 0,
    height: 110,
    selection: ko.observable(),
    clicktimer: undefined,
    collapse: ko.observable(false),
    
    
    activate: function() {
      var self = this;
      self.dataset(ko.utils.unwrapObservable(require('viewmodels/visualisation').dataset));
    
      app.on('allele:selected').then(function(d) {
        if (d) {
          self.scrollToOffset(d);
        }
      });
    },
    
    attached: function() {
      var self = this;

      self.render();
    },
    
    render: function() {
      var self = this;
      

      var ref = this.dataset().reference_sequence();

      // Cut sample string in blocks
      var sample_sequence = this.dataset().sample_sequence();
      var reference_sequence = this.dataset().reference_sequence();
      var last_end_position = 1;
      self.segments = [];
      self.dataset().allele_description().forEach(function(mutation) {
        // make the string between last known end position and start of new mutation if > 0
        var unmodified_length = mutation.reference_start() - last_end_position;
        
        // If there is a piece of unmodified text between the new and the last mutation, add this to the array
        if (unmodified_length > 0) {
          self.segments.push({
            reference_start: last_end_position,
            reference_end: mutation.reference_start()-1,
            annotations_top: null,
            annotations_bottom: null,
            type: null,
            unmodified_sequence: reference_sequence.substr(last_end_position-1, unmodified_length),
            old_pieces: [],
            new_pieces: [],
          });
        }
        
        var mutation_sample_length = mutation.sample_end() - mutation.sample_start();
        self.segments.push({
          reference_start: mutation.reference_start(),
          reference_end: mutation.reference_end(),
          annotations_top: null,
          annotations_bottom: null,
          type: mutation.type(),
          unmodified_sequence: null,
          old_pieces: mutation.pieces().filter(function(segment) { return segment.type() == 'del'; }),
          new_pieces: mutation.pieces().filter(function(segment) { return segment.type() != 'del'; }),
          hgvs: mutation.hgvs(),
          mutation: mutation,
        });
        
        // Due to differences in position numbering for different mutation types, this statement is needed
        if (mutation.type() == 'ins' || mutation.type() == 'dup') {
          last_end_position = mutation.reference_end();
        } else {
          last_end_position = mutation.reference_end();
        }
      });
      
      if (last_end_position < reference_sequence.length) {
        self.segments.push({
          reference_start: last_end_position,
          reference_end: reference_sequence.length-1,
          annotations_top: null,
          annotations_bottom: null,
          type: null,
          unmodified_sequence: reference_sequence.substr(last_end_position-1),
          old_pieces: [],
          new_pieces: [],
        });
        
      }
      
      var d3_segments = d3.select('.graph-string').classed('animation-state', false).selectAll('.string-segment').data(self.segments);
      var new_segments = d3_segments.enter().append('div')
        .attr('class', function(d) { return 'string-segment '+ d.type; })
        .on('click', self.handleClick.bind(self));
      
      var annotations = new_segments.append('div').filter(function(d) { return d.type === null;}).attr('class', 'annotation-top');
      annotations.append('span').attr('class', 'start').text(function(d) { return d.reference_start; });
      annotations.append('span').attr('class', 'end').text(function(d) { return d.reference_end; });
      
      annotations = new_segments.append('div').filter(function(d) { return d.type !== null;}).attr('class', 'annotation-hgvs').text(function(d) { return d.hgvs; });
      
      var old_pieces = new_segments.append('div').filter(function(d) { return d.type !== null;})
        .attr('class', 'old')
        .selectAll('.piece')
        .data(function(d) { return d.old_pieces; })
        .enter()
          .append('div')
          .attr('class', function(d) { return 'piece piece-'+d.type()+(d.reverse() ? ' inverse': ''); });
      
      old_pieces.append('div').attr('class', 'piece-annotation').text(function(d) { return (d.type() == 'trans') ? 'transposition from ' + d.origin_start() + ' – '+ d.origin_end() : ''; });
      old_pieces.append('div').attr('class', 'piece-str').text(function(d) { return (0 && d.sequence().length > 50) ? d.sequence().substr(0, 4) + '⋯' + d.sequence().substr(-4,4) : d.sequence(); });

      new_segments.append('div').filter(function(d) { return d.type === null;})
        .attr('class', 'unmodified').text(function(d) { return ( d.unmodified_sequence.length > 1500) ? d.unmodified_sequence.substr(0, 4) + '⋯' + d.unmodified_sequence.substr(-4,4) : d.unmodified_sequence; });

      var new_pieces = new_segments.append('div').filter(function(d) { return d.type !== null;})
        .attr('class', 'new')
        .selectAll('.piece')
        .data(function(d) { return d.new_pieces; })
        .enter().append('div')
            .attr('class', function(d) { return 'piece piece-'+d.type()+(d.reverse() ? ' inverse': ''); });

      new_pieces.append('div').attr('class', 'piece-str').text(function(d) { return (d.sequence().length > 50) ? d.sequence().substr(0, 4) + '⋯' + d.sequence().substr(-4,4) : d.sequence(); });
      new_pieces.append('div').filter(function(d) { return (d.type() == 'trans'); }).attr('class', 'piece-annotation').text(function(d) { return 'transposition from ' + d.origin_start() + ' – '+ d.origin_end(); });
      
      d3.select('.graph-string').append('span').text('');
    },
    
    /**
     * Scroll the view to a certain offset
     * @method scrollToOffset
     * @param {Mutation} d the mutation to scroll into view
     */
    scrollToOffset: function(d) {
      var self = this;
      console.log('self segments', self.segments);
      var domNode = d3.select('.graph-string').selectAll('.string-segment').data(self.segments).filter(function(item) { return item.hgvs == ko.utils.unwrapObservable(d.hgvs); }).node();
      // domNode.scrollIntoView();
      $('.graph-string').animate( {scrollLeft: $(domNode).offset().left - 120 + $('.graph-string').scrollLeft()}, {duration: 350, easing: 'swing'});
    },
    
    /**
     * Action handler for click on event.
     * Fires a new global event that is received by all viewmodels
     * 
     * @param {Mutation} d The mutation
     * @method handleClick
     */
    handleClick: function(d) {
      console.log('click', d);
      if (d.mutation) app.trigger('allele:selected', d.mutation);
    },
    
  };
});



