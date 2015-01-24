/**
 * Generate paths from the allele description.
 * 
 * @module layouts/sequence
 * @exports d3.layout.sequence
 */
define(['d3', 'knockout'], function (d3, ko) {  
  if ( typeof(d3.layout.sequence) !== 'undefined')
    return d3.layout.sequence;
    
  /** 
   * Constructor of new sequence layouter.
   * 
   * @constructor
   */
  d3.layout.sequence = function() {
    
    var rx, sx, ry, sy, growthSize = 5;
    
    /**
     * Generate paths from allele description.
     * **At the moment Knockout observables are assumed**.
     * 
     * @method sequence()
     * @param {Mixed} data Data to convert containing knockout observables.
     * @param {function} [index] callback returning key for data
     */
    function sequence(data, index) {
      // Convert data to array of groups
      data.forEach(function(mutation,i) {
        sequence.createPlaceholder(mutation);
        mutation.pieces().forEach(function(segment, i2) {
          sequence.createShape(mutation, segment);
        });
      });
      return data;
    }
    
    sequence.createPlaceholder = function(mutation) {
      mutation.shape = sequence.makePolygon(
        ko.utils.unwrapObservable(mutation.reference_start),
        ko.utils.unwrapObservable(mutation.reference_end),
        ko.utils.unwrapObservable(mutation.sample_start),
        ko.utils.unwrapObservable(mutation.sample_end),
        ko.utils.unwrapObservable(false)
      );
    };

    /**
     * Set / get accessor of pixel x coordinate on reference axis. 
     * 
     * @method sequence.rx
     * @param   {function}  [x] Callback returning pixel x coordinate on reference value for data.
     * @return  {function|object}    Returns callback function if `x` was empty, returns self otherwise.
     */
    sequence.rx = function(x) {
      if (!arguments.length) return rx;
      rx = x;
      return sequence;
    };
    
    /**
     * Set / get accessor of pixel y coordinate of reference axis. 
     * 
     * @method sequence.ry
     * @param   {double}  [y] Pixel y coordinate of reference axis
     * @return  {double|object}    Returns callback function if `y` was empty, returns self otherwise.
     */
    sequence.ry = function(y) {
      if (!arguments.length) return ry;
      ry = y;
      return sequence;
    };
    
    /**
     * Set / get accessor of pixel x coordinate on sample axis. 
     * 
     * @method sequence.sx
     * @param   {function}  [x] Callback returning pixel x coordinate on sample value for data.
     * @return  {function|object}    Returns callback function if `x` was empty, returns self otherwise.
     */
    sequence.sx = function(x) {
      if (!arguments.length) return sx;
      sx = x;
      return sequence;
    };
    
    /**
     * Set / get accessor of pixel y coordinate of sample axis. 
     * 
     * @method sequence.sy
     * @param   {double}  [x] Pixel y coordinate of sample axis
     * @return  {double|object}    Returns callback function if `x` was empty, returns self otherwise.
     */
    sequence.sy = function(x) {
      if (!arguments.length) return sy;
      sy = x;
      return sequence;
    };
    
    /**
     * Set / get accessor of pixel x coordinate on sample axis. 
     * 
     * @method sequence.hoverGrowth
     * @param   {double}  [s]   Size that each object should be extended for an hover class.
     * @return  {double|object}    Returns growth size if `s` was empty, returns self otherwise.
     */
    sequence.hoverGrowth = function(s) {
      if (!arguments.length) return hoverGrowth;
      hoverGrowth = s;
      return sequence;
    };
    
    /** 
     * Calculate the intersection between the lines (r1 to s1) and (r2 to s2)
     *
     * @method sequence.calculateIntersection
     * @param  {Double}  r1  position on reference-axis in **pixel coordinates** of first diagonal.
     * @param  {Double}  s1  position on sample-axis in **pixel coordinates** of first diagonal.
     * @param  {Double}  r1  position on reference-axis in **pixel coordinates** of first diagonal.
     * @param  {Double}  s1  position on sample-axis in **pixel coordinates** of first diagonal.
     * @return  {Object}  x and y pixel coordinates of point of intersection.
     */
    sequence.calculateIntersection = function(r1, s1, r2, s2) {
      // Assume the default `yt` and `yr` for the calculation of `dy`.
      var dy = Math.abs(sy - ry);
      
      // Write first line as `y1 = a * x + c`
      var a = dy / (r1 - s1);
      var c = - a * s1;
      
      // Write second line as `y2 = b * x + d`
      var b = dy / (r2 - s2);
      var d = - b * s2;
      
      // Calculate the intersection between the two lines.
      var intersection_x = (d - c) / (a - b);
      var intersection_y = sy - (a * ( intersection_x ) + c);
      
      return {x: intersection_x, y: intersection_y};
    };
    
    /**
     * Create a SVG Path description out of two ranges on the reference and sample axis
     * 
     * @method sequence.makePolygon
     * @param   {Double}  r1  First value on reference axis
     * @param   {Double}  r2  Second value on reference axis
     * @param   {Double}  s1  First value on sample axis
     * @param   {Double}  s2  Second value on sample axis
     * @param   {Boolean}  reverse  Set to true if the polygon should be reversed, results in hourglass shape
     * @return  {String}    SVG 'd' path description of polygon
     */
    sequence.makePolygon =  function(r1, r2, s1, s2, reverse) {
      // if r2 / s2 are missing, it's operating on a single position. Set range to 1.
      if (!r2) r2 = r1 + 1;
      if (!s2) s2 = s2 + 1;
      
      // Shift all elements one half to the left, such that the axis labels will indicate the center of a letter instead of the space between.
      // r1 += -0.5;
      // r2 += -0.5;
      // s1 += -0.5;
      // s2 += -0.5;
      
      // Determine the order of the values and convert them to pixel coordinates.
      r_left = rx((r1 < r2) ? r1 : r2);
      r_right = rx((r1 < r2) ? r2 : r1);
      t_left = sx((s1 < s2) ? s1 : s2);
      t_right = sx((s1 < s2) ? s2 : s1);
      
      if (reverse) {
        return ' M ' + r_left + ' ' + ry + 
               ' L ' + r_right + ' ' + ry + 
               ' L ' + t_left + ' ' + sy + 
               ' L ' + t_right + ' ' + sy + 
               ' Z ';

      } else {
        return ' M ' + r_left + ' ' + ry +
               ' L ' + r_right + ' ' + ry +
               ' L ' + t_right + ' ' + sy +
               ' L ' + t_left + ' ' + sy +
               ' Z ';
      }
    };
    
    /**
     * Create a SVG Path description out of two ranges on the reference and sample axis, 
     * growing the region in the x-direction in a fixed distance perpendicular to each line.
     * 
     * @method sequence.makeExtendedPolygon
     * @param   {Double}  r1  First value on reference axis
     * @param   {Double}  r2  Second value on reference axis
     * @param   {Double}  s1  First value on sample axis
     * @param   {Double}  s2  Second value on sample axis
     * @param   {Boolean}  reverse  Set to true if the polygon should be reversed, results in hourglass shape
     * @return  {String}    SVG 'd' path description of polygon
     */
    sequence.makeExtendedPolygon = function(r1, r2, s1, s2, reverse) {
      if (typeof(reverse) === 'undefined') reverse = false;
      var ext = growthSize; // size of the extension
      var ext_left;
      var ext_right;
      var intersection_left;
      var intersection_right;
      var path;
      
      // If `r2` and `s2` are missing, it's operating on a single position. Set range to 1.
      if (!r2) r2 = r1 + 1;
      if (!s2) s2 = s2 + 1;
      
      // Shift all elements one half to the left, such that the axis labels will indicate the center of a letter instead of the space between.
      // r1 += -0.5;
      // r2 += -0.5;
      // s1 += -0.5;
      // s2 += -0.5;
      
      // Determine the order of the variables (r1, r2) and (s1,s2) and convert to pixel-coordinates.
      r_left = rx((r1 < r2) ? r1 : r2);
      r_right = rx((r1 < r2) ? r2 : r1);
      t_left = sx((s1 < s2) ? s1 : s2);
      t_right = sx((s1 < s2) ? s2 : s1);
      
      // If the polygon is reversed, create a shape like an hourglass.
      if (reverse) {
        ext_left = ext / Math.sin( Math.atan2( (sy - ry),  (t_right - r_left )));
        ext_right = ext / Math.sin( Math.atan2( (sy - ry),  (t_left - r_right )));
        
        // Calculate two intersections, one shifted by `ext_left` of the hourglass and one shifted by `ext_right` of the hourglass.
        intersection_left  = sequence.calculateIntersection(r_left - ext_left, t_right - ext_left, r_right - ext_right, t_left - ext_right);
        intersection_right = sequence.calculateIntersection(r_left + ext_left, t_right + ext_left, r_right + ext_right, t_left + ext_right);
        
        path = ' M ' + (r_left - ext_left) + ' ' + ry +
                   ' L ' + (r_right + ext_right) + ' ' + ry +
                   ' L ' + intersection_right.x + ' ' + intersection_right.y +
                   ' L ' + (t_right + ext_left) + ' ' + sy +
                   ' L ' + (t_left - ext_right) + ' ' + sy +
                   ' L ' + intersection_left.x + ' ' + intersection_left.y +
                   ' Z ';
        return path;
      } else {
        // If the polygon is not reversed, just output a normal polygon.
        ext_left = ext / Math.sin( Math.atan2( (sy - ry),  (t_left - r_left )));
        ext_right = ext / Math.sin( Math.atan2( (sy - ry),  (t_right - r_right )));
        
        // polygon is not inverted, so draw normal polygon:
        return ' M ' + (r_left - ext_left) + ' ' + ry +
               ' L ' + (r_right + ext_right) + ' ' + ry +
               ' L ' + (t_right + ext_right) + ' ' + sy +
               ' L ' + (t_left - ext_left) + ' ' + sy +
               ' Z ';
      }
    };
    
    /**
     * Create a path object from data
     *
     * @method createShape
     * @param {Object} d Data for path
     * @return {Object} Object containing path-information
     */
    sequence.createShape = function(mutation, segment) {
      if (segment.trans()) {
        segment.shape = sequence.makePolygon(
          ko.utils.unwrapObservable(segment.origin_start),
          ko.utils.unwrapObservable(segment.origin_end),
          ko.utils.unwrapObservable(segment.sample_start),
          ko.utils.unwrapObservable(segment.sample_end),
          ko.utils.unwrapObservable(segment.reverse)
        );
        
        segment.hover = sequence.makeExtendedPolygon(
          ko.utils.unwrapObservable(segment.origin_start),
          ko.utils.unwrapObservable(segment.origin_end),
          ko.utils.unwrapObservable(segment.sample_start),
          ko.utils.unwrapObservable(segment.sample_end),
          ko.utils.unwrapObservable(segment.reverse)
        );
        
      } else if (segment.type() == 'ins' && (mutation.type() == 'ins' || mutation.type() == 'delins')) {
        segment.shape = sequence.makePolygon(
          ko.utils.unwrapObservable(mutation.reference_end),
          ko.utils.unwrapObservable(mutation.reference_end),
          ko.utils.unwrapObservable(segment.sample_start),
          ko.utils.unwrapObservable(segment.sample_end),
          ko.utils.unwrapObservable(segment.reverse)
        );
        
        segment.hover = sequence.makeExtendedPolygon(
          ko.utils.unwrapObservable(mutation.reference_end),
          ko.utils.unwrapObservable(mutation.reference_end),
          ko.utils.unwrapObservable(segment.sample_start),
          ko.utils.unwrapObservable(segment.sample_end),
          ko.utils.unwrapObservable(segment.reverse)
        );
        
      } else if (segment.type() == 'del' && (mutation.type() == 'delins')) {
        segment.shape = sequence.makePolygon(
          ko.utils.unwrapObservable(mutation.reference_start),
          ko.utils.unwrapObservable(mutation.reference_end),
          ko.utils.unwrapObservable(segment.sample_start),
          ko.utils.unwrapObservable(segment.sample_start),
          ko.utils.unwrapObservable(segment.reverse)
        );
        
        segment.hover = sequence.makeExtendedPolygon(
          ko.utils.unwrapObservable(mutation.reference_start),
          ko.utils.unwrapObservable(mutation.reference_end),
          ko.utils.unwrapObservable(segment.sample_start),
          ko.utils.unwrapObservable(segment.sample_start),
          ko.utils.unwrapObservable(segment.reverse)
        );
      } else if (segment.type() == 'del' && (mutation.type() == 'inv')) {
        // Don't show the deleted part of an invertion.
      } else {
        segment.shape = sequence.makePolygon(
          ko.utils.unwrapObservable(mutation.reference_start),
          ko.utils.unwrapObservable(mutation.reference_end),
          ko.utils.unwrapObservable(segment.sample_start),
          ko.utils.unwrapObservable(segment.sample_end),
          ko.utils.unwrapObservable(segment.reverse)
        );
        
        segment.hover = sequence.makeExtendedPolygon(
          ko.utils.unwrapObservable(mutation.reference_start),
          ko.utils.unwrapObservable(mutation.reference_end),
          ko.utils.unwrapObservable(segment.sample_start),
          ko.utils.unwrapObservable(segment.sample_end),
          ko.utils.unwrapObservable(segment.reverse)
        );
      }
            
      return segment;
    };  
    
    return sequence;
  };
  
  return d3.layout.sequence;
});
