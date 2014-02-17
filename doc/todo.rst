Todo list
=========

These are some general todo notes. More specific notes can be found by
grepping the source code for ``Todo``.

.. seealso::

   `Mutalyzer Trac -- Active tickets <https://humgenprojects.lumc.nl/trac/mutalyzer/report/2>`_
     Users can file tickets on the Mutalyzer Trac website.

   `Mutalyzer GitLab -- Open issues <https://git.lumc.nl/mutalyzer/mutalyzer/issues>`_
     Some issues are recorded in the Mutalyzer GitLab project.

- Improve the web interface design :)
- Test all uses of mkstemp().
- Use naming conventions for modules Crossmap, Db, File, GenRecord, Retriever
  and Scheduler.
- Use standard logging module, with rotating functionality. Race conditions
  on the log file are probably a problem in the current setup.
  Instead of that rotating, we could also use logrotate:
  http://serverfault.com/questions/55610/logrotate-and-open-files
- Setup continuous integration. Currently, I'm most impressed with Hudson.
  - http://hudson-ci.org/
  - http://www.rhonabwy.com/wp/2009/11/04/setting-up-a-python-ci-server-with-hudson/
  Or perhaps Jenkins.
  - http://jenkins-ci.org/
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
