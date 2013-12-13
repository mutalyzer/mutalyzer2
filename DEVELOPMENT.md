Mutalyzer development
=====================

Development of Mutalyzer happens on the GitLab server:
https://git.lumc.nl/mutalyzer/mutalyzer


Coding style
------------

In general, try to follow the [PEP 8](http://www.python.org/dev/peps/pep-0008)
guidelines for Python code and
[PEP 257](http://www.python.org/dev/peps/pep-0257/) for docstrings.


Unit tests
----------

The unit tests depend on a running batch daemon, webserver, and SOAP
web service:

    sudo /etc/init.d/mutalyzer-batchd start
    sudo /etc/init.d/apache2 start

Now run the tests with:

    nosetests -v

Or, if you are in a hurry, skip the long-running tests with:

    MUTALYZER_QUICK_TEST=1 nosetests -v


Working with feature branches
-----------------------------

New features are best implemented in their own branches, isolating the work
from unrelated developments. In fact, it's good practice to **never work
directly on the master branch** but always in a separate branch. For this
reason, the master branch on the GitLab server is locked. Feature branches can
be merged back into master via a **merge request** in GitLab.

Before starting work on your feature, create a branch for it:

    git checkout -b your-feature

Commit changes on this branch. If you're happy with it, push to GitLab:

    git push origin your-feature -u

Now create a merge request to discuss the implementation with your
colleagues. This might involve adding additional commits which are included in
the merge request by pushing your branch again:

    git commit
    git push

You may also be asked to rebase your branch on the master branch if it has
changed since you started your work. This will require a forced push:

    git fetch
    git rebase origin/master
    git push -f

If the work is done, a developer can merge your branch and close the merge
request. After the branch was merged you can safely delete it:

    git branch -d your-feature


Release management
------------------

The current Mutalyzer version is recorded in `mutalyzer/__init__.py`. See the
comments in that file for more info of the versioning scheme.

On the event of a new release, the following is done:

    emacs mutalyzer/__init__.py

Update `__date__`, remove `dev` from `__version_info__` and set `RELEASE` to
`True`.

    git commit -am 'Bump version to 2.0.beta-XX'
    git tag -a 'mutalyzer-2.0.beta-XX'
    git push --tags

    emacs mutalyzer/__init__.py

Set `__version_info__` to a new version ending with `'dev'` and set `RELEASE`
to `FALSE`.

    git commit -am 'Open development for 2.0.beta-YY'

Be sure to upgrade your installations to the new version as described in the
INSTALL file (e.g. `sudo python setup.py develop` for development checkouts).


Development notes
-----------------

Todo list:

- Improve the web interface design :)
- Test all uses of mkstemp().
- Use naming conventions for modules Crossmap, Db, File, GenRecord, Retriever
  and Scheduler.
- Use standard logging module, with rotating functionality. Race conditions
  on the log file are probably a problem in the current setup.
  Instead of that rotating, we could also use logrotate:
  http://serverfault.com/questions/55610/logrotate-and-open-files
- Setup continuous integration. Currently, I'm most impressed with Hudson.
  http://hudson-ci.org/
  http://www.rhonabwy.com/wp/2009/11/04/setting-up-a-python-ci-server-with-hudson/
  Or perhaps Jenkins.
  http://jenkins-ci.org/
- Migrate Javascript to JQuery.
- I think in the long run, the Output object is not really the way to go. It
  obscures the control flow. The logging part should use the standard logging
  module. The data gathering by the Output object is probably better handled
  by explicitely returning data objects from functions.
- Migrate from TAL to a more mondern and maintained Python template library,
  for example jinja.
- Develop a large test suite.
- Create a web interface url to watch the progress of a batch job.
- Create web services for the batch jobs (steal ideas from Jeroen's DVD
  web service).
- Use virtualenv?
- Use SQLAlchemy?
- Password for MySQL user.
- In deployment, remove old versions of Mutalyzer package?
- Check for os.path.join vulnerabilities.
- Use a standard solution for the database migrations in extras/migrations.
- Use something like Sphinx to generate development documentation from code.
- There are some problems with the batch architecture, especially that there
  cannot be multiple workers without synchronisation problems.
  Good read: http://news.ycombinator.com/item?id=3002861
  Suggestion: http://celeryproject.org/
- Have a normal 404 page.
- Maintenance (and/or read-only) mode.
- Cleanup this document.
- Be more explicit in all the type of descriptions we don't currently support.

Code style guide:

- Follow PEP 8 (code) and PEP 257 (docstrings).
    http://www.python.org/dev/peps/pep-0008/
    http://www.python.org/dev/peps/pep-0257/
  Read the Google Python Style guide:
    http://google-styleguide.googlecode.com/svn/trunk/pyguide.html
- Use Epydoc style documentation in docstrings.
- End class and method definitions with their name as comment.
- Executables are in the bin/ directory.
- For examples, check established Python projects:
    http://code.djangoproject.com/browser/django/trunk
    http://twistedmatrix.com/trac/browser/trunk
    https://github.com/webpy/webpy
    https://github.com/mitsuhiko/jinja2
    https://bitbucket.org/mramm/tg-21/src
    http://bazaar.launchpad.net/~bzr-pqm/bzr/bzr.dev/files
    https://github.com/ask/celery
- A lot of code does not yet adhere to these points, this is an ongoing
  effort.

Obsoleted features:

- On eu.liacs.nl:
  /etc/apache2/mods-enabled/rewrite.load contains a rewrite rule that converts
  "Variant_info.php" to "Variant_info".
  When all LOVD versions are above 2.0-23, this rule can be deleted and the
  rewrite module can be disabled.
- In the Variant_info() function a substitution on error messages is
  performed.
  When all LOVD versions are above 2.0-23, this check can be deleted.


Dependencies
------------

Mutalyzer depends on the following (Debian/Ubuntu) packages:

- mysql-server     >= 5.1
- python           >= 2.6
- python-mysqldb   >= 1.2.2
- python-biopython >= 1.54
- python-pyparsing >= 1.5.0
- python-configobj >= 4.4.0
- python-magic     >= 5.04-2
- python-psutil    >= 0.1.3-1
- python-xlrd      >= 0.6.1-2
- python-daemon    >= 1.5.5
- python-soappy    >= 0.12.0-2
- python-suds      >= 0.3.9-1

The web and SOAP interfaces depend on the following packages:

- apache2             >= 2.2.11
- libapache2-mod-wsgi >= 2.8
- python-webpy        >= 0.33
- python-rpclib       >= 2.8.0-beta
- python-simpletal    >= 4.1-6

Automatic remote deployment depends on Fabric:

- fabric >= 0.9.0-2

The unit tests depend on the following packages:

- python-nose    >= 0.11
- python-webtest >= 1.2.3

As of 2011-08-23, snakefood reports the following imports from the Mutalyzer
source code (excluding the standard library imports):

    Bio
    MySQLdb
    SOAPpy
    configobj
    daemon
    fabric
    lockfile
    lxml
    magic
    nose
    pyparsing
    setuptools
    simpletal
    rpclib
    suds
    web
    webtest
    xlrd
