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

To run the unit tests with `pytest`_, just run ``py.test``.


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

A normal version number takes the form X.Y.Z where X is the major version, Y
is the minor version, and Z is the patch version. Development versions take
the form X.Y.Z.dev where X.Y.Z is the closest future release version.

Note that this scheme is not 100% compatible with `SemVer`_ which would
require X.Y.Z-dev instead of X.Y.Z.dev but `compatibility with setuptools
<http://peak.telecommunity.com/DevCenter/setuptools#specifying-your-project-s-version>`_
is more important for us. Other than that, version semantics are as described
by SemVer.

Releases are available from the GitLab git repository as tags. Additionally,
the latest release is available from the `release` branch.

.. note:: Older Mutalyzer version numbers took the form 2.0.beta-X where X was
   incremented on release.


Release procedure
^^^^^^^^^^^^^^^^^

Releasing a new version is done as follows:

1. Make sure the section in the ``CHANGES`` file for this release is
   complete and there are no uncommitted changes.

   .. note::

    Commits since release X.Y.Z can be listed with ``git log vX.Y.Z..`` for
    quick inspection.

2. Update the ``CHANGES`` file to state the current date for this release
   and edit ``mutalyzer/__init__.py`` by updating `__date__` and removing the
   ``dev`` value from `__version_info__`.

   Commit and tag the version update::

       git commit -am 'Bump version to X.Y.Z'
       git tag -a 'vX.Y.Z'

3. Push to the GitLab repository (assuming the remote name is `gitlab` and you
   are working on the `master` branch::

       git push gitlab master
       git push gitlab master:release --tags

4. Add a new entry at the top of the ``CHANGES`` file like this::

       Version X.Y.Z+1
       ---------------

       Release date to be decided.

   Increment the patch version and add a ``dev`` value to `__version_info__`
   in ``mutalyzer/__init__.py`` and commit these changes::

       git commit -am 'Open development for X.Y.Z+1'


.. _pytest: http://pytest.org/
.. _PEP 8: http://www.python.org/dev/peps/pep-0008/
.. _PEP 257: http://www.python.org/dev/peps/pep-0257/
.. _SemVer: http://semver.org/
