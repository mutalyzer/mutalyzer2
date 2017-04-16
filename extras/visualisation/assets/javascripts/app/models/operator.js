/**
 * Operator model
 * 
 * @module models/operator
 * @exports Operator
 */
define(["knockout"], function (ko) {
	var Operator = function Operator(data) {
		self = this;
		self.name = ko.observable(data.name);
		self.class = ko.observable(data.class);
		self.enabled = ko.observable(true);
		self.type = ko.observable();
		self.realType = ko.computed(function(){
			return self.type();
		});
	};
	return Operator;
});