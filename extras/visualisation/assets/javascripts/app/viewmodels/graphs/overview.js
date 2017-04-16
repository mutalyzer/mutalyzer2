/**
 * Overview viewmodel.
 * 
 * @module viewmodels/graphs/overview
 */
define(['d3', 'durandal/app', 'knockout'], function (d3, app, ko) {
  
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
    height: 150,
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
      app.on('allele:hovered').then(self.displayOperationInfo.bind(self));
    },
    
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
      
      self.displayOperationInfo(d);
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
      console.log(d);
      self.svg.classed('has-hover', d !== undefined);
         
      self.svg.selectAll(".allele.hover").classed('hover', false);
      
      if (d) {
        self.svg.selectAll(".allele").data([d], function(d) { return d.hgvs(); }).classed('hover', true);
      }
      
      app.trigger('allele:hovered', d);
    },
    
    /**
     * Action handler for doubleclick event
     * Seperate single and double click events
     * 
     * @param {Mutation} d The mutation
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
      self.svg.classed('has-hover', false);
      self.svg.classed('has-active', false);
      

      // Doubleclick is on already selected element, undo zoom.
      if (d === undefined || !force && this.selection() == d) {
        self.svg.select('.brush-reference').call(self.brushReference.extent([0,0]));
        self.svg.select('.brush-sample').call(self.brushSample.extent([0,0]));
        self.selection(undefined);
        self.highlightMutation(undefined);
        return;
      }

      // Grow ranges 5% on each axis.
      var ext_reference = Math.max(10, Math.min(self.reference.domain()[1]*0.05, 0.5 * (1 + d.r2() - d.r1())));
      var ext_sample = Math.max(10, Math.min(self.sample.domain()[1]*0.05, 0.5 * (1 + d.s2() - d.s1())));

      self.svg.select('.brush-reference').call(self.brushReference.extent([
          Math.max(d.r1() - ext_reference, 0), 
          Math.min(d.r2() + ext_reference, self.reference.domain()[1])
        ]));
        
      self.svg.select('.brush-sample').call(self.brushSample.extent([
          Math.max(d.s1()-ext_sample, 0),
          Math.min(d.s2() + ext_sample, self.sample.domain()[1])
        ]));

      self.selection(d);
      self.highlightMutation(d);
      
      return false;
    },
    
    /**
     * Called by Durandal upon attachment of this view to the DOM
     * @method attached
     */
    attached: function() {
      var self = this;
      self.render();
      console.log('resizing zoom2');
      
      app.on('window:resize', function() {
        self.render();
        self.selectionChanged(self.selection(), true);
      });
      
      self.dataset().allele_description.subscribe(function(dataset) {
        self.refreshData();
      });      
    },
    
    /**
     * Display a hovering box with extra information.
     * @param {Mutation} d Display information on this mutation.
     * @method displayOperationInfo
     */
    displayOperationInfo: function(d) {
      var self = this;
      if (!d) {
        self.tooltip.transition().duration(500).style('opacity', 1e-6);
        return;
      }
      
      self.tooltip.transition().duration(200)
        .style('opacity', 1)
        .attr('class', 'graph-tooltip ' + d.type())
        .style('left', ( (self.sample(d.s1()) + self.sample(d.s2())) / 2) + 'px');
      self.tooltip
        .html('<strong>'+d.type()+'</strong>'+d.hgvs())
        .style('bottom', '-20px');
    },
    
    /**
     * Callback method for brushes.
     * @method brushReferenceStart
     */
    brushReferenceStart: function() {
      var self = this;
      self.svg.classed("selecting", true);
    },
    
    /**
     * Callback method for brushes.
     * @method brushReferenceMove
     */
    brushReferenceMove: function() {
      var self = this;
      if (self.brushReference.extent()[0] !== self.brushReference.extent()[1]) {
        app.trigger('graph:domain:selection:change', { domain: self.brushReference.extent(), axis: 'reference', sender: self});
      }
    },
    
    /**
     * Callback method for brushes.
     * @method brushReferenceEnd
     */
    brushReferenceEnd: function() {
      var self = this;
      if (self.brushReference.extent()[0] == self.brushReference.extent()[1]) {
        self.brushReference.extent(self.reference.domain());
      }
      app.trigger('graph:domain:selection:change', { domain: self.brushReference.extent(), axis: 'reference', sender: self});
      self.svg.classed("selecting", !d3.event.target.empty());
    },
    
    /**
     * Callback method for brushes.
     * @method brushSampleStart
     */
    brushSampleStart: function() {
      var self = this;
      self.svg.classed("selecting", true);
    },
    
    /**
     * Callback method for brushes.
     * @method brushSampleMove
     */
    brushSampleMove: function() {
      var self = this;
      if (self.brushSample.extent()[0] !== self.brushSample.extent()[1]) {
        app.trigger('graph:domain:selection:change', { domain: self.brushSample.extent(), axis: 'sample', sender: self});
      }
    },
    
    /**
     * Callback method for brushes.
     * @method brushSampleEnd
     */
    brushSampleEnd: function() {
      var self = this;
      if (self.brushSample.extent()[0] == self.brushSample.extent()[1]) {
        self.brushSample.extent(self.sample.domain());
      }
      app.trigger('graph:domain:selection:change', { domain: self.brushSample.extent(), axis: 'sample', sender: self});
      self.svg.classed("selecting", !d3.event.target.empty());
    },
    
    /**
     * Render the static parts of the visualisation.
     * @method render
     */
    render: function() {
      var self = this;
      self.svg = d3.select('.graph-overview svg');
      $('.graph-overview svg').empty();
      self.width = $('.graph-overview').width();
      
      self.margin = {top: 40, right: 20, bottom: 30, left: 60};
      
      self.reference = d3.scale.linear().range([self.margin.left, self.width - self.margin.right]).domain([1, self.dataset().reference_sequence().length]);
      self.sample   = d3.scale.linear().range([self.margin.left, self.width - self.margin.right]).domain([1, self.dataset().sample_sequence().length]);

      self.yr = self.margin.top + 1;
      self.yt = self.height - self.margin.bottom - 1;
      
      self.layout = new d3.layout.sequence();
      self.layout.rx(self.reference).sx(self.sample).ry(self.yr).sy(self.yt);
      
      self.rAxis = d3.svg.axis().scale(self.reference).orient('top').ticks(20).tickFormat(d3.format('d'));
      self.tAxis = d3.svg.axis().scale(self.sample).orient('bottom').ticks(20).tickFormat(d3.format('d'));
                  
      self.svg.attr('width', this.width).attr('height', this.height);
      
      self.svg.append("g")
          .attr("class", "axis t")
          .attr("transform", "translate(0," + (self.height - self.margin.bottom)+ ")")
          .call(self.tAxis)
        .append("text")
          .attr('x', self.margin.left-5).attr('y', 0).attr('dy', '0.5em')
          .style("text-anchor", "end")
          .text("Sample");
          
      self.svg.append("g")
          .attr("class", "axis r")
          .attr("transform", "translate(0," + self.margin.top+ ")")
          .call(self.rAxis)
        .append("text")
          .attr('x', self.margin.left-5).attr('y', 0).attr('dy', '0.5em')
          .style("text-anchor", "end")
          .text("Reference");
          
      self.brushReference = d3.svg.brush()
        .x(self.reference)
        .on("brushstart", self.brushReferenceStart.bind(self))
        .on("brush", self.brushReferenceMove.bind(self))
        .on("brushend", self.brushReferenceEnd.bind(self));
        
      self.brushSample = d3.svg.brush()
        .x(self.sample)
        .on("brushstart", self.brushSampleStart.bind(self))
        .on("brush", self.brushSampleMove.bind(self))
        .on("brushend", self.brushSampleEnd.bind(self));
      
        
      self.arc = d3.svg.arc()
        .outerRadius(10)
        .startAngle(0)
        .endAngle(function(d, i) { return i ? -Math.PI : Math.PI; });

      self.brushReferenceG = self.svg.append("g")
        .attr("class", "brush brush-reference")
        .call(self.brushReference);

      self.brushReferenceG.selectAll(".resize").append("path")
        .attr("transform", "translate(0," +  (self.margin.top-10)+ ")")
        .attr("d", self.arc);

      self.brushReferenceG.selectAll("rect")
        .attr("transform", "translate(0," +  (self.margin.top-20)+ ")")        
        .attr("height", 20);
        
        
      self.brushSampleG = self.svg.append("g")
        .attr("class", "brush brush-sample")
        .call(self.brushSample);

      self.brushSampleG.selectAll(".resize").append("path")
        .attr("transform", "translate(0," +  (self.height - self.margin.bottom + 10)+ ")")
        .attr("d", self.arc);

      self.brushSampleG.selectAll("rect")
        .attr("transform", "translate(0," +  (self.height - self.margin.bottom )+ ")")        
        .attr("height", 20);

      self.tooltip = d3.select(".graph-overview .graph-tooltip").style('opacity', 0);
      
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
      
      var groups = self.svg.selectAll(".allele").data(data, function(d) { return d.hgvs(); });
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
          .append('path').attr('class', 'placeholder').attr('d', function(d) { return (ko.utils.unwrapObservable(d.hasTrans)) ? ko.utils.unwrapObservable(d.shape) : ''; });
          
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
