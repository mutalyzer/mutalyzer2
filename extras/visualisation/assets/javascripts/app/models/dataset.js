/**
 * Dataset model
 * 
 * @module models/dataset
 * @exports Dataset
 */
define(["knockout"], function (ko) {
	var Dataset = function Dataset(data) {
		self = this;
		self.id = ko.observable(data.id);
		self.name = ko.observable(data.name);
		self.date = ko.observable(data.date);
		self.url = ko.observable(data.url);
		self.length = ko.observable(data.length);
		self.description = ko.observable(data.description);
		
		self.reference_sequence = ko.observable("");
		self.target_sequence = ko.observable("");
		self.allele_description = ko.observableArray([]);
	};
	return Dataset;
});