## Mutalyzer DNA Visualisation

### Usage
#### Installation

Prerequisites: 

The project depends on a lot of dependencies. The most important are:
- [NodeJS](https://github.com/joyent/node/wiki/Installing-Node.js-via-package-manager) 
- [Mimosa](http://mimosa.io) is used as a build tool;
- [Durandal](http://durandaljs.com) is a Single Page Application framework;
- [D3](http://d3js.org) is used to render the visualisations;
- [Groc](http://nevir.github.io/groc/) is used to generate documentation.


Install [NodeJS](https://github.com/joyent/node/wiki/Installing-Node.js-via-package-manager) first. After that, run the following commands to install

```bash
npm install -g bower mimosa groc
npm install
bower install
```

#### Starting development server
The setup contains a simple development server. The makefile runs the Mimosa tool. This tool watches changes in the code and will refresh the browser automatically. The project is served from the `public` folder.

```bash
make start
```

#### Building the project
Building the project creates a minimised and optimised version of the code. All dependencies will be concatenated in one javascript file. 

```bash
make dist
```

#### Creating the documentation

The code is documented with docblock tags. The documentation can be converted with the following command:

```bash
groc assets/javascripts/app/**/*.js README.md
```


### Project layout

    assets/
      javascripts/
        app/                   # Root of real application
          backend/backend.js   # Interaction with backend
          layouts/sequence.js  # Creates shapes from data
          models/              # Models for representing data
          viewmodels/          # Contains js files with view logic
          views/               # HTML part of views
        vendor/                # Third party javascript
     stylesheets/              # styling of the app
       vendor/                 # Third party CSS / LESS
    bower_components/          # third party frontend libraries
    dist/                      # Compiled version of the app
    doc/                       # Default target dir for Groc
    node_modules/              # NPM server side libraries
    routes/index.js            # Used for running local server
    views/index.hbs            # Used for running local server


#### Fixed dataset
The dataset is fixed in the `activate` function in  `visualisation.js`. Also, the routes in `shell.js` must be changed.

-- 
Copyright 2015 Landscape - http://www.wearelandscape.nl