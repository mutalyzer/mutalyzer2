/**
 * Mutation model
 * 
 * @module models/mutation
 * @exports Mutation
 */
define(["knockout", "komapping"], function (ko, komapping) {
	var Mutation = function Mutation(data) {
		var self = this;
		self.type = ko.observable(data.type);
		self.reverse = ko.observable(data.reverse);
		self.pieces = komapping.fromJS(data.pieces);
    
    // Add a reference back to this object for all pieces
    self.pieces().forEach(function(piece) { piece.operator = self; });
    
    // In some datasets the hgvs was missing. Create another indicator for the list
		self.hgvs = ko.observable(data.hgvs || data.type+'/'+data.reference_start+'/'+data.reference_end+'/'+data.sample_start+'/'+data.sample_end);
		
		self.reference_start = ko.observable(data.reference_start);
		self.reference_end = ko.observable(data.reference_end);
		self.sample_start = ko.observable(data.sample_start);
		self.sample_end = ko.observable(data.sample_end);
		self.reverse = ko.observable(data.reverse);
		self.nested = ko.observable(data.nested);
		
		self.hasTrans = ko.computed(function() {
			return self.pieces().some(function(piece) { return piece.trans(); });
		});
		
		self.r1 = ko.computed(function() {
			return Math.min.apply(Math, self.pieces().map(function(d) { return d.origin_start(); }).filter(function(d) { return d > 0; }).concat([self.reference_start()]));
		});
		
		self.r2 = ko.computed(function() {
			return Math.max.apply(Math, self.pieces().map(function(d) { return d.origin_end(); }).concat([self.reference_end()]));
		});
		
		self.s1 = ko.computed(function() {
			return self.sample_start();
		});
		
		self.s2 = ko.computed(function() {
			return self.sample_end();
		});
	};
  
	return Mutation;
});