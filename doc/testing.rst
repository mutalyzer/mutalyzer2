.. highlight:: none

.. _testing:

Testing
=======


Unit tests
----------

We use `pytest`_ for the unit tests. To run them, just type ``py.test`` from
the Mutalyzer source directory.

.. note:: The Mutalyzer package must be installed before running the unit
          tests.

Tests are `run automatically on Travis CI
<https://travis-ci.org/LUMC/mutalyzer>`_ for each pull request and push on
GitHub.


Testing the web services
------------------------

To ease testing the web services during development, some simple web service
client scripts are included in the Mutalyzer source tree::

    $ cd extras/soap-tools
    $ ./info.py
    Version: 2.0.5
    Version parts: 2, 0, 5
    Release date: 16 Dec 2014
    Nomenclature version: 2.0
    Nomenclature version parts: 2, 0
    Server name: res-muta-app01
    Contact e-mail: humgen@lumc.nl
    $

They simply call one of the web service functions and print the result. You
may have to change the server location defined at the top of these scripts.

.. note:: One of the scripts, ``run_batch_job.py``, provides an easy way to
          run a batch job from the command line. Some additional notes are
          available for `running this on a Windows machine
          <https://gist.github.com/jfjlaros/482fe9f0397e554ed29f>`_.


.. _pytest: http://pytest.org/
