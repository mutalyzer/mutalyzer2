define(['durandal/app', 'knockout', 'backend/backend'], function (app, ko, backend) {
	return {
		displayName: 'Datasets',
		datasets: ko.observableArray([]),
		selectDataset: function(dataset) {
			app.showMessage('Loading dataset '+dataset.name());
		},
    
    /**
     * Called by Durandal upon activating the ViewController
     * @method activate
     */
		activate: function () {
			var self = this;
			backend.getDatasets().then(function(results) {
				self.datasets(results);
			}, function(err) { console.error(err); });
		},
		attached: function() {
		}
	};
});