define(['plugins/http', 'jquery', 'models/operator', 'models/dataset', 'models/mutation',  'komapping', 'knockout'], function (http, $, Operator, Dataset, Mutation, komapping, ko) {
  return {
    _datasets: [],
  
    /**
     * Convert mutation input from json-format to internal format
     * 
     * @method newMapMutations
     * @param {Array} allele Original list of mutations
     * @param {String} reference Reference DNA-string
     * @param {String} sample Sample DNA-string
     * @return {Array} List of hidden mutations
     */
    mapMutations: function(allele, reference, sample)
    {
      var self = this;
      var offset = 0;

      // Loop over all mutations 
      return allele.map(function(mutation,i) {
        
        // For readability of debug output, delete these variables;
        delete mutation.start_aa;
        delete mutation.end_aa;
        delete mutation.start_offset;
        delete mutation.end_offset;
        delete mutation.shift;
        delete mutation.term;
        delete mutation.sample_start_offset;
        delete mutation.sample_end_offset;
        delete mutation.dna;
        delete mutation.hgvs_length;
        
        // Convert all ranges to array-based indexing
        // For insertions, convert `end = start+1` to `end = start`
        if (mutation.type == 'ins' || mutation.type == 'dup') {
          mutation.end = mutation.start;
          mutation.sample_end++;
        } else if (mutation.type == 'del') {
          mutation.sample_start = mutation.sample_end;
          mutation.end++;
        } else {
          mutation.sample_end++;
          mutation.end++;
        }
        
        // There is an error in the input file format for the `inv` mutation type. Correct this while the file format is wrong
        if (mutation.type == 'inv') {
          mutation.inserted[0].reverse = true;
          mutation.inserted[0].sequence = reference.substr(mutation.sample_start-1, mutation.sample_end - mutation.sample_start);
        }
        
        // If the mutation has a deleted part
        if (~['del', 'inv', 'delins', 'subst'].indexOf(mutation.type)) {
          mutation.deleted[0].start = mutation.start;
          mutation.deleted[0].end = mutation.end;
          // All pieces should have their own `sample_start` and `sample_end` parameters;
          mutation.deleted[0].sample_start = mutation.sample_start;
          mutation.deleted[0].sample_end = mutation.sample_start;
        } else {
          // The mutation type doesn't support deleted pieces so delete them.
          mutation.deleted = [];
        }
        
        // If the mutation has an inserted part
        if (~['ins', 'delins', 'subst', 'inv', 'dup'].indexOf(mutation.type)) {
          offset = 0;

          // Loop over all inserted pieces.
          mutation.inserted = mutation.inserted.map(function(piece) {
            // Not only correct the `end` notation to index-based, also insert the `sample_start` and `sample_end`
            // parameters. These are used to
            if (piece.type == 'trans') {
              piece.end++;
              piece.sample_start = mutation.sample_start + offset;
              offset += piece.end - piece.start;
              piece.sample_end = mutation.sample_start + offset;
              
              // Fetch sequence from sample string, because it can be inverted
              piece.sequence = reference.substr(piece.sample_start, piece.sample_end - piece.sample_start);
              
            } else {
              piece.sample_start = mutation.sample_start + offset;
              offset += piece.sequence.length;
              piece.sample_end = mutation.sample_start + offset;
            }
            
            // In the original json-format, the `sample_end` of the mutation has the wrong value for mutations with
            // multiple segments. While this is not fixed, this should be corrected.
            mutation.sample_end = piece.sample_end;
            return piece;
          });
        } else {
          // The mutation type doesn't support inserted pieces so delete them.
          mutation.inserted = [];
        }
        
        offset = 0;
        var mutation_out =  {
          type: mutation.type,
          hgvs: mutation.hgvs,
          sample_start: mutation.sample_start,
          sample_end: mutation.sample_end,
          reference_start: mutation.start,
          reference_end: mutation.end,
          nested: (mutation.inserted.length > 1),
          pieces:
            mutation.deleted.map(function(piece, i2) {
              return {
                type: 'del',
                mutationtype: mutation.type,
                trans: false,
                reverse: false,
                sequence: piece.sequence, 
                origin_start: null,
                origin_end: null,
                reference_start: mutation.start,
                reference_end: mutation.end,
                // Deletions result in something between two elements
                sample_start: piece.sample_start,
                sample_end: piece.sample_start
              };
            })
          .concat(
            mutation.inserted.map(function(piece,i2) {
              var segment_out =  {
                type: piece.type,
                mutationtype: mutation.type,
                trans: (piece.type == 'trans'),
                reverse: piece.reverse,
                sequence: piece.sequence || '--ERROR--',
                origin_start: piece.start,
                origin_end: piece.end,
                reference_start: mutation.start,
                reference_end: mutation.start, // new elements are inserted at beginning
                sample_start: piece.sample_start,
                sample_end: piece.sample_end,
              };
              
              return segment_out;
            })
          ),
        };
        
        return mutation_out;
      });    
    },
    
    /**
     * Retrieve a dataset from the server
     * 
     * @method getDataset
     * @param {Integer} id of the dataset in the datasets.json
     * @return {Promise} Promise that resolves when dataset is loaded and mapped.
     */
    getDataset: function(id)
    {
      var self = this;
      var dfd = $.Deferred();

      self.getDatasets().then(function(datasets) {
        var dataset = ko.utils.arrayFirst(datasets, function(item) { return (item.id() == id); });
        http.get(dataset.url()).then(function(results) {
          
          // As we now have the correct format of the JSON input, the mapping is not needed anymore.
          results.allele_description = self.mapMutations(results.allele_description, 'r'+results.reference_sequence, 's'+results.sample_sequence);
          var r = komapping.fromJS(results, {
            allele_description: {
              create: function(options) {
                return new Mutation(options.data);
              }
            }
          }, dataset);
          dfd.resolve(r);
        }, function(err) {
          console.error('There was an error loading the dataset', err);
          dfd.reject(err);
        });
      });
      return dfd.promise();
    },
    
    /**
     * Retrieve list of datasets from the server.
     * 
     * @method getDatasets
     * @return {Promise} A promise that resolves when the list of datasets is loaded
     */
    getDatasets: function()
    {
      var self = this;
      var dfd = $.Deferred();
      
      if (self._datasets.length === 0) {
        return http.get('datasets.json').then(function(results) {
          var r = $.map(results, function(item) { return new Dataset(item); });
          dfd.resolve(r);
          return dfd.promise();
        });
      } else {
        dfd.resolve(self._datasets);
        return dfd.promise();
      }
    },
    
    /**
     * Retrieve list of operators from the server.
     * 
     * @method getOperators
     * @return {Promise} A promise that resolves when the list of operators is loaded
     */
    getOperators: function()
    {
      return http.get('operators.json').then(function(results)
      {
        var dfd = $.Deferred();
        var r = $.map(results, function(item) { return new Operator(item); });
        dfd.resolve(r);
        
        return dfd.promise();
      });
    }
  };
});
