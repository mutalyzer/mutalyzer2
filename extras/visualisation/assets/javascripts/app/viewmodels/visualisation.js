define(['durandal/app', 'knockout', 'komapping', 'backend/backend'], function (app, ko, komapping, backend) {
  return {
    displayName: 'Visualisation',
    operators: ko.observableArray([]),
    dataset: ko.observable(),
    selectedOperation: ko.observable(),
    referenceDomain: { begin: ko.observable(), end: ko.observable() },
    sampleDomain: { begin: ko.observable(), end: ko.observable() },
    inAnimation: ko.observable(false),
    
    /**
     * Called by Durandal upon activating the ViewController
     * @method activate
     * @return {Promise} A promise that is resolved after the dataset is loaded.
     */
    activate: function (id) {
      var promise
      var self = this;
      
      // use this to fix the dataset at 4
      if (!id) id = 4;
      
      self.selectOperation.bind(self, self);

      // Register trigger on selection change event
      app.on('allele:selected').then(self.selectionChanged.bind(self));
      
      promise = backend.getDataset(id).then(
        function(results) {
          self.dataset(results);
          self.referenceDomain.begin(1);
          self.referenceDomain.end(results.reference_sequence().length);
          self.sampleDomain.begin(1);
          self.sampleDomain.end(results.sample_sequence().length);
        },
        function(err) {
          app.showMessage('Error loading the dataset', 'Error');
        }
      );
      
      backend.getOperators().then(function(results) {
        komapping.fromJS(results, {}, self.operators);
      }, function(err) { console.error(err); });
    
      // Register trigger on zoomed range changes.
      app.on('graph:domain:selection:change', function(ev) {
        if (ev.axis == 'reference') {
           self.referenceDomain.begin(Math.floor(ev.domain[0]));
           self.referenceDomain.end(Math.ceil(ev.domain[1]));
        } else {
           self.sampleDomain.begin(Math.floor(ev.domain[0]));
           self.sampleDomain.end(Math.ceil(ev.domain[1]));
        }
      });

      // Register trigger for resetting zoom.
      app.on('graph:domain:selection:reset', function(ev) {
        self.referenceDomain.begin(1);
        self.referenceDomain.end(self.dataset().reference_sequence().length);
        self.sampleDomain.begin(1);
        self.sampleDomain.end(self.dataset().sample_sequence().length);
      });
      
      return promise;
    },
    
    /**
     * Called by Durandal upon attaching the view to the DOM
     * @method attached
     */
    attached: function() {
      //
      $(window).on('resize', function() {
        app.trigger('window:resize');
      }); 
      
    },
      
    /**
     * Animate the mutations
     * @method playAlleleDescription
     */
    playAlleleDescription: function() {
      var self = this;
      var string;
      var allele;
      var promise_start;
      var promise;
      
      self.inAnimation(true);
      
      string = require('viewmodels/stringanimation');
      string.resetString();

      self.resetZoom();
      self.selectedOperation(null);
      allele = self.dataset().allele_description();
      self.dataset().allele_description([]);      

      promise_start = $.Deferred();
      promise = allele.reduce(function(promise, d) {
        return promise.then(function() {
          self.dataset().allele_description.push(d);
        }).then(function(result) {
          app.trigger('movie:step', d);
          self.selectedOperation(d);
          return d.pieces().reduce(function(promise2,piece) {
            return promise2.then(function(result) {
              return string.animate(piece);
            });
          }, promise);
        });
      }, promise_start);
      
      promise.then(function(){
        string.resetString();
        // string.inAnimation(false);
        self.inAnimation(false);
        self.dataset().allele_description(allele);
        //   setTimeout(self.playAlleleDescription.bind(self), 2000);
      });
      
      promise_start.resolve();
    },
    
    /**
     * Called when changing the selection
     * This method is linked to the `allele:selected` event and responds to 
     * changes from within the view as well as from other views.
     * 
     * @param {Mutation} d The mutation
     * @method selectionChanged
     */
    selectionChanged: function(d) {
      var self = this;
      self.selectedOperation((self.selectedOperation() == d) ? undefined : d);
    },
    
    /**
     * Action handler for clicking on a mutation in the sidebar.
     * Fires a new global event that is received by all viewmodels
     * 
     * @param {Mutation} d The mutation
     * @method selectOperation
     */
    selectOperation: function(d) {
      var self = this;
      app.trigger('allele:selected', d);
    },
    
    /**
     * Trigger handler for reset zoom event
     * 
     * @method resetZoom
     */
    resetZoom: function() { 
      var self = this;
      app.trigger('graph:domain:selection:reset', {axis: 'reference', sender: self});
      app.trigger('graph:domain:selection:reset', {axis: 'sample', sender: self});
    }
  };
});