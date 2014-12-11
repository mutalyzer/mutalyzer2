.. highlight:: none

.. _contributing:

Contributing
============

Contributions to Mutalyzer are very welcome! They can be feature requests, bug
reports, bug fixes, unit tests, documentation updates, or anything els you may
come up with.

Development of Mutalyzer happens on GitHub:
https://github.com/LUMC/mutalyzer


Coding style
------------

In general, try to follow the `PEP 8`_ guidelines for Python code and `PEP
257`_ for docstrings.


Installation
------------

As a developer, you probably want to install the Mutalyzer package in
development mode. This will allow you to edit files directly in the source
directory without having to reinstall.

Please refer to :ref:`install` for general installation instructions. For
development mode installation, instead of using ``python setup.py install``,
use::

    python setup.py develop


Creating a pull request
-----------------------

Contributions are most welcome as GitHub pull requests. If you're familiar
with the typical GitHub pull request workflow, you can skip this section.

New features are best implemented in their own branches, isolating the work
from unrelated developments. In fact, it's good practice to *never work
directly on the master branch* but always in a separate branch. When your work
is ready, a feature branch can be merged back into master via a *pull request*
in GitHub.

Before starting your work, fork the Mutalyzer repository to your own namespace
in GitHub and work from this fork. Before starting work on your feature,
create a branch for it::

    git clone https://github.com/<you>/mutalyzer.git && cd mutalyzer
    git checkout -b your-feature

Commit changes on this branch. If you're happy with it, push to GitHub::

    git push origin your-feature -u

Now create a pull request to discuss the implementation with the Mutalyzer
maintainers. This might involve adding additional commits which are included
in the pull request by pushing your branch again::

    git commit
    git push

If the work is done, a maintainer can merge your branch and close the pull
request. After the branch was merged you can safely delete it::

    git branch -d your-feature


.. _PEP 8: http://www.python.org/dev/peps/pep-0008/
.. _PEP 257: http://www.python.org/dev/peps/pep-0257/
