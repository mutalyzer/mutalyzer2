.. highlight:: none

.. _development:

Development
===========

Development of Mutalyzer happens on GitLab:
https://git.lumc.nl/mutalyzer/mutalyzer


Contributing
------------

Contributions to Mutalyzer are very welcome! They can be feature requests, bug
reports, bug fixes, unit tests, documentation updates, or anything els you may
come up with.


Coding style
------------

In general, try to follow the `PEP 8`_ guidelines for Python code and `PEP
257`_ for docstrings.


Unit tests
----------

To run the unit tests with `nose`_, just run ``nosetests -v``.


Working with feature branches
-----------------------------

New features are best implemented in their own branches, isolating the work
from unrelated developments. In fact, it's good practice to *never work
directly on the master branch* but always in a separate branch. For this
reason, the master branch on the GitLab server is locked. Feature branches can
be merged back into master via a *merge request* in GitLab.

Before starting work on your feature, create a branch for it::

    git checkout -b your-feature

Commit changes on this branch. If you're happy with it, push to GitLab::

    git push origin your-feature -u

Now create a merge request to discuss the implementation with your
colleagues. This might involve adding additional commits which are included in
the merge request by pushing your branch again::

    git commit
    git push

You may also be asked to rebase your branch on the master branch if it has
changed since you started your work. This will require a forced push::

    git fetch
    git rebase origin/master
    git push -f

If the work is done, a developer can merge your branch and close the merge
request. After the branch was merged you can safely delete it::

    git branch -d your-feature


Versioning
----------

All version numbers for recent Mutalyzer releases take the form 2.0.beta-X
where X is incremented on release. Pre-release (or development) version
numbers take the form 2.0.beta-X.dev where 2.0.beta-X is the closest future
release version.

Note that we are planning a switch to `SemVer`_.

.. A normal version number takes the form X.Y.Z where X is the major version, Y
   is the minor version, and Z is the patch version. Development versions take
   the form X.Y.Z.dev where X.Y.Z is the closest future release version.

   Note that this scheme is not 100% compatible with `SemVer`_ which would
   require X.Y.Z-dev instead of X.Y.Z.dev but `compatibility with setuptools
   <http://peak.telecommunity.com/DevCenter/setuptools#specifying-your-project-s-version>`_
   is more important for us. Other than that, version semantics are as described
   by SemVer.

   Releases are `published at PyPI <https://pypi.python.org/pypi/wiggelen>`_ and
   available from the GitHub git repository as tags.


Release procedure
^^^^^^^^^^^^^^^^^

Releasing a new version is done as follows:

1. Make sure the section in the ``CHANGES`` file for this release is
   complete and there are no uncommitted changes.

   .. note::

    Commits since release 2.0.beta-X can be listed with ``git log
    mutalyzer-2.0.beta-X..`` for quick inspection.

2. Update the ``CHANGES`` file to state the current date for this release
   and edit ``mutalyzer/__init__.py`` by updating `__date__`, removing the
   ``dev`` value from `__version_info__` and setting `RELEASE` to `True`.

   Commit and tag the version update::

       git commit -am 'Bump version to 2.0.beta-X'
       git tag -a 'mutalyzer-2.0.beta-X'
       git push --tags

3. Add a new entry at the top of the ``CHANGES`` file like this::

       Version 2.0.beta-Y
       ------------------

       Release date to be decided.

   Set `__version_info__` to a new version ending with ``dev`` and set
   `RELEASE` to `True` in ``mutalyzer/__init__.py``. Commit these changes::

       git commit -am 'Open development for 2.0.beta-Y'


.. _nose: https://nose.readthedocs.org/
.. _PEP 8: http://www.python.org/dev/peps/pep-0008/
.. _PEP 257: http://www.python.org/dev/peps/pep-0257/
.. _SemVer: http://semver.org/
