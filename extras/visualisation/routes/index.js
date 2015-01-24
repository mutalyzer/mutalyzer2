var index = function(config) {
  cachebust = ''
  if (process.env.NODE_ENV !== "production") {
    cachebust = "?b=" + (new Date()).getTime()
  }

  var options = {
    reload: config.liveReload.enabled,
    optimize: config.isOptimize != null ? config.isOptimize : false,
    cachebust: cachebust
  };

  return function(req, res) {
    res.render("index", options);
  };
};

exports.index = index;