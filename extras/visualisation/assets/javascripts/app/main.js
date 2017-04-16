requirejs.config({
	paths: {
		'text': '../vendor/requirejs-text/text',
		'knockout': '../vendor/knockout.js/knockout',
		'jquery': '../vendor/jquery/jquery',
		'bootstrap': '../vendor/bootstrap/bootstrap',
		'durandal':'../vendor/durandal',
		'plugins' : '../vendor/durandal/plugins',
		'transitions' : '../vendor/durandal/transitions',
    'd3' : '../vendor/d3/d3',
		'backend' : 'backend',
		'komapping': "../vendor/knockout-mapping/knockout.mapping",
		'layouts': 'layouts',
	},
	shim: {
		'bootstrap': {
				deps: ['jquery'],
				exports: 'jQuery'
		}
	}
});

define(function(require) {
	var app = require('durandal/app'),
		viewLocator = require('durandal/viewLocator'),
		system = require('durandal/system');

	//>>excludeStart("build", true);
	system.debug(true);

	//>>excludeEnd("build");

	app.title = 'Mutalyzer DNA Visualisation';

	app.configurePlugins({
		router: true,
		dialog: true,
		widget: true
	});

	app.start().then(function() {
		//Replace 'viewmodels' in the moduleId with 'views' to locate the view.
		//Look for partial views in a 'views' folder in the root.
		viewLocator.useConvention();

		//Show the app by setting the root view model for our application with a transition.
		app.setRoot('viewmodels/shell');
	});
});