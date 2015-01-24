define(['d3', 'durandal/app', 'knockout', 'layouts/sequence'], function (d3, app, ko) {

  /**
   * Used to reorder the elements in the visualisation
   * @method d3.selection.prototype.moveToFront
   */  
  d3.selection.prototype.moveToFront = function() {
    return this.each(function(){
      this.parentNode.appendChild(this);
    });
  };
    
  return {
    root: null,
    dataset: ko.observable(),
    width: 0,
    height: 200,
    selection: ko.observable(),
    clicktimer: undefined,
    
    /**
     * Called by Durandal upon activating the ViewController
     * @method activate
     */
    activate: function() {
      var self = this;
      self.dataset(ko.utils.unwrapObservable(require('viewmodels/visualisation').dataset));

      // Register event handlers for selection changes
      app.on('allele:selected').then(self.selectionChanged.bind(self));

      // Register handler for hovering a mutation
      app.on('allele:hovered').then(self.hoverMutation.bind(self));
      
      // Register event handler for changes in the domain.
      app.on('graph:domain:selection:change', function(ev) {
        if (ev.axis == 'reference') {
           self.reference.domain(ev.domain);
          self.rAxis.ticks(Math.min(30, ev.domain[1] - ev.domain[0]));
           self.svg.select('.axis.r').call(self.rAxis).call(function(selection) {
            selection.selectAll('.tick text').attr('transform', 'translate(' + (self.reference(1.5) - self.reference(1)) + ', 0)');
          });
           self.refreshData();
        } else if (ev.axis == 'sample') {
           self.sample.domain(ev.domain);
          self.tAxis.ticks(Math.min(30, ev.domain[1] - ev.domain[0]));
           self.svg.select('.axis.t').call(self.tAxis).call(function(selection) {
            selection.selectAll('.tick text')
            .attr('transform', 'translate(' + (self.sample(1.5) - self.sample(1)) + ', 0)');
          });
           self.refreshData();
        }
      });

      // Register event handler for a reset of the selection.
      app.on('graph:domain:selection:reset', function(ev) {
        if (ev.axis == 'reference') {
          self.reference.domain(self.fullReference.domain());
          self.rAxis.ticks(Math.min(20, self.fullReference.domain()[1] - self.fullReference.domain()[0]));
           self.svg.select('.axis.r').call(self.rAxis);
           self.refreshData();
        } else {
          self.sample.domain(self.fullTarget.domain());
          self.tAxis.ticks(Math.min(20, self.fullTarget.domain()[1] - self.fullTarget.domain()[0]));
           self.svg.select('.axis.t').call(self.tAxis);
           self.refreshData();
        }
      });
      
      // Register event handler for animation
      app.on('movie:step', function(d) {
        self.selectionChanged(d);
      });
    },
  
    /**
     * Action handler for doubleclick event
     * Seperate single and double click events
     * 
     * @param {Mutation} d 
     * @method handleClick
     */
    handleClick: function(d) {
      app.trigger('allele:selected', d);
    },
    
    /**
     * Called when changing the selection
     * This method is linked to the `allele:selected` event and responds to 
     * changes from within the view as well as from other views.
     * 
     * @param {Mutation} d The mutation
     * @method selectionChanged
     */
    selectionChanged: function(d, force) {
      var self = this;
      
      // doubleclick is on already selected element, undo zoom.
      if (!force && this.selection() == d || d === undefined) {
        self.reference.domain(self.fullReference.domain());
        self.rAxis.ticks(Math.min(20, self.fullReference.domain()[1] - self.fullReference.domain()[0]));
        self.svg.select('.axis.r').call(self.rAxis);
        self.sample.domain(self.fullTarget.domain());
        self.tAxis.ticks(Math.min(20, self.fullTarget.domain()[1] - self.fullTarget.domain()[0]));
        self.svg.select('.axis.t').call(self.tAxis);
        self.selection(undefined);
        self.refreshData();
        self.highlightMutation(undefined);
        return;
      }
      
      self.selection(d);

      // Grow ranges 5% on each axis. Use a minimal zoomlevel:
      var ext_reference = Math.max(10, Math.min(self.fullReference.domain()[1]*0.05, 0.5 * (1 + d.r2() - d.r1())));
      var ext_sample = Math.max(10, Math.min(self.fullTarget.domain()[1]*0.05, 0.5 * (1 + d.s2() - d.s1())));
      
      var ref_domain = [
        Math.max(d.r1() - ext_reference, 0), 
        Math.min(d.r2() + ext_reference, self.fullReference.domain()[1])
      ];
      
      var sample_domain = [
        Math.max(d.s1()-ext_sample, 0),
        Math.min(d.s2() + ext_sample, self.fullTarget.domain()[1])
      ];
      
      self.reference.domain(ref_domain);
      self.sample.domain(sample_domain);
      
      self.rAxis.ticks(Math.min(30, ref_domain[1] - ref_domain[0]));
      self.tAxis.ticks(Math.min(30, sample_domain[1] - sample_domain[0]));
      
      self.svg.select('.axis.r').call(self.rAxis).call(function(selection) {
        selection.selectAll('.tick text').attr('transform', 'translate(' + (self.reference(1.5) - self.reference(1)) + ', 0)');
      });
      
      self.svg.select('.axis.t').call(self.tAxis).call(function(selection) {
        selection.selectAll('.tick text').attr('transform', 'translate(' + (self.sample(1.5) - self.sample(1)) + ', 0)');
      });
      
      self.selection(d);
      self.refreshData();
      self.highlightMutation(d);
      
      return false;
    },
  
    /**
     * Action handler for click on event.
     * Fires a new global event that is received by all viewmodels
     * 
     * @param {Mutation} d The mutation
     * @method selectMutation
     */
    
    /**
     * Highlight a mutation
     * This method is linked to the `allele:hovered` event and responds to 
     * changes from within the view as well as from other views.
     * 
     * @param {Mutation} d The mutation
     * @method highlightMutation
     */
    highlightMutation: function(d) {
      var self = this;
         
      self.svg.classed('has-active', d !== undefined);
      self.svg.selectAll(".allele.active").classed('active', false);
      
      if (d) {
        self.svg.selectAll(".allele").data([d], function(d) { return d.hgvs(); }).classed('active', true);
      }

    },
    
    /**
     * Action handler for hover event.
     * Fires a new global event that is received by all viewmodels
     * 
     * @param {Mutation} d The mutation
     * @method hover
     */
    hoverMutation: function(d) {
      var self = this;
      self.svg.classed('has-hover', d !== undefined);
         
      self.svg.selectAll(".allele.hover").classed('hover', false);
      
      if (d) {
        self.svg.selectAll(".allele").data([d], function(d) { return d.hgvs(); }).classed('hover', true);
      }
    },
    
    /**
     * Action handler for hover event.
     * Fires a new global event that is received by all viewmodels
     * 
     * @param {Mutation} d The mutation
     * @method hover
     */
    hover: function(d) {
      var self = this;
      app.trigger('allele:hovered', d);
    },
    
    /**
     * Called by Durandal upon attachment of this view to the DOM
     * @method attached
     */
    attached: function() {
      var self = this;
      self.render();

      app.on('window:resize', function() {
        self.render();
        self.selectionChanged(self.selection(), true);
      });
      self.dataset().allele_description.subscribe(function(dataset) {
        self.refreshData();
      });
    },
    
    /**
     * Render the static parts of the visualisation.
     * @method render
     */
    render: function() {
      var self = this;
      self.svg = d3.select('.graph-zoom svg');
      $('.graph-zoom svg').empty();
      self.width = $('.graph-zoom').width();
      
      self.margin = {top: 40, right: 20, bottom: 30, left: 60};
      
      self.reference = d3.scale.linear().range([0, self.width - self.margin.right - self.margin.left]).domain([1, self.dataset().reference_sequence().length]);
      self.sample   = d3.scale.linear().range([0, self.width - self.margin.right - self.margin.left]).domain([1, self.dataset().sample_sequence().length]);
      self.fullReference = d3.scale.linear().range([0, self.width - self.margin.right - self.margin.left]).domain([1, self.dataset().reference_sequence().length]);
      self.fullTarget   = d3.scale.linear().range([0, self.width - self.margin.right - self.margin.left]).domain([1, self.dataset().sample_sequence().length]);

      self.yr = self.margin.top + 1;
      self.yt = self.height - self.margin.bottom - 1;
      
      self.layout = new d3.layout.sequence();
      self.layout.rx(self.reference).sx(self.sample).ry(self.yr).sy(self.yt);
      
      self.rAxis = d3.svg.axis().scale(self.reference).orient('top').ticks(20).outerTickSize(0).tickFormat(d3.format('d'));
      self.tAxis = d3.svg.axis().scale(self.sample).orient('bottom').ticks(20).outerTickSize(0).tickFormat(d3.format('d'));
                  
      self.svg.attr('width', this.width).attr('height', this.height);
      self.svg.append('g').attr("transform", "translate(" + (self.margin.left)+ ", 0)").append('svg').attr('class', 'mutations').attr('width', self.width - self.margin.left - self.margin.right);

      
      self.svg.append("g")
          .attr("class", "axis t")
          .attr("transform", "translate("+self.margin.left+"," + (self.height - self.margin.bottom)+ ")")
          .call(self.tAxis)
        .append("text")
          .attr('x', -5).attr('y', 0).attr('dy', '0.5em')
          .style("text-anchor", "end")
          .text("Sample");
      
      
      self.svg.append("g")
          .attr("class", "axis r")
          .attr("transform", "translate("+self.margin.left+"," + self.margin.top+ ")")
          .call(self.rAxis)
          .call(function(selection) {
            selection.selectAll('.major text')
            .attr('transform', 'translate(' + self.reference(0.5) + ', 0)');
          })
        .append("text")
          .attr('x', -5).attr('y', 0).attr('dy', '0.5em')
          .style("text-anchor", "end")
          .text("Reference");
          
      self.refreshData();
    },
    
    /**
     * Render dynamic parts of the data.
     * This function is called each time the data or parameters are changed.
     * @method refreshData
     */
    refreshData: function() {
      var self = this;
      var data = self.layout(self.dataset().allele_description());
      var groups = self.svg.select('.mutations').selectAll(".allele").data(data, function(d) { return d.hgvs(); });
      groups.enter()
        .append('g')
          .attr("class", function(d) { return 'allele '+ko.utils.unwrapObservable(d.type); })
          .on('click', self.handleClick.bind(self))
          .on("mouseover", function() {
            var sel = d3.select(this);
            sel.moveToFront();
            sel.each(function(d) { self.hover(d); });
          })
          .on('mouseout', function() { self.hover(undefined); })
          .append('path').attr('class', 'placeholder');
      
          groups.select('.placeholder').attr('d', function(d) { return (ko.utils.unwrapObservable(d.hasTrans)) ? ko.utils.unwrapObservable(d.shape) : ''; });
      
      var shapes = groups.selectAll('g').data(function(d) { return d.pieces(); });
      var newshapes = shapes.enter().append('g').attr('class', 'shape');
      newshapes.append('path').attr('class', 'hover-area').attr('d', function(d) { return ko.utils.unwrapObservable(d.hover); });
      newshapes.append('path').attr('class', function(d) { return 'shape '+ko.utils.unwrapObservable(d.type); });
      
      shapes.select('.hover-area').attr('d', function(d) { return ko.utils.unwrapObservable(d.hover); });
      shapes.select('.shape').attr('d', function(d) { return ko.utils.unwrapObservable(d.shape); });
      
      groups.exit().remove();
    },
    
  };
});
