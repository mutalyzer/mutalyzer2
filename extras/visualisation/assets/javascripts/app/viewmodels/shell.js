define(['plugins/router'], function (router) {
	return {
		router: router,
    
    /**
     * Called by Durandal upon activating the ViewController
     * In this place the routes are defined. 
     * @method activate
     */
		activate: function () {
			router.map([
          // Change this route to for example 'datasets' in order to support dynamic dataset loading.
					{ route: '', title:'Visualisation', moduleId: 'viewmodels/visualisation', nav: true},
					{ route: 'datasets', title:'Datasets', moduleId: 'viewmodels/datasets', nav: true},
					{ route: 'visualisation/:id', title:'Visualisation', moduleId: 'viewmodels/visualisation', hash: '#visualisation', nav: true},
				]).buildNavigationModel();

			return router.activate();
		}
	};
});