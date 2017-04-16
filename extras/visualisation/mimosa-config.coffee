exports.config =
  
  minMimosaVersion:'2.0.0'

  modules: [
    'server'
    'require'
    'minify-js'
    'minify-css'
    'live-reload'
    'less'
    'combine'
    'requirebuild-include'
    'requirebuild-textplugin-include'
    'bower'
    'web-package'
    'csslint'
    'jshint'
    'copy'
  ]

  watch:
    javascriptDir: 'javascripts/app'

  requireBuildTextPluginInclude:
    pluginPath: 'text'
    extensions: ['html']

  requireBuildInclude:
    folder:"javascripts"
    patterns: ['app/**/*.js', 'vendor/durandal/**/*.js']

  bower:
    bowerDir:
      clean: true
      path: ".mimosa/bower_components"
    copy:
      enabled: true
      trackChanges: true
      strategy: 'packageRoot'
      mainOverrides:
        # "knockout.js":["knockout.js","knockout-2.3.0.debug.js"]
        "knockout-mapping": [
          "knockout.mapping.js"
        ]
        "bootstrap": [
          "dist/css/bootstrap.css"
          "dist/js/bootstrap.js"
        ]
        "font-awesome": [
          { font: "../../font" }
          "css/font-awesome.css"
          "css/font-awesome-ie7.css"
        ]
        "durandal": [
          {
            img: "../../images"
            js: "durandal"
            css: "durandal"
          }
        ]
        "requirejs":["require.js"]
        "requirejs-text":["text.js"]

  combine:
    folders: [
      {
        folder:'stylesheets'
        output:'stylesheets/styles.css'
        order: [
          'vendor/bootstrap/bootstrap.css'
          'vendor/font-awesome/font-awesome.css'
          'vendor/durandal/durandal.css'
          'starterkit.css'
        ]
      }
    ]

  server:
    defaultServer:
      enabled: true
      onePager: true
    port: 3001
    views:
      compileWith: 'handlebars'
      extension: 'hbs'

  require:
    verify:
      enabled: true
    optimize:
      overrides:
        paths:
          d3 : '../vendor/d3/d3'
          backend : 'backend'
          komapping: "../vendor/knockout-mapping/knockout.mapping"
          layouts: 'layouts'
        name: '../vendor/almond-custom'
        inlineText: true
        stubModules: ['text']
        pragmas:
          build: true

