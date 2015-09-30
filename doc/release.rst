.. highlight:: none

.. _releases:

Releases
========


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

Releases are available from the GitHub git repository as tags. Additionally,
the latest release is available from the `release` branch.

.. note:: Older Mutalyzer version numbers took the form 2.0.beta-X where X was
   incremented on release.


Release procedure
-----------------

Releasing a new version is done as follows. This assumes remote `github` is
the upstream Mutalyzer repository and you have push rights there.

1. Start a release branch and make sure the section in the ``CHANGES.rst``
   file for this release is complete::

       git checkout -b release-X.Y.Z
       git add CHANGES.rst
       git commit -m 'Update changelog'

   .. note::

    Commits since release X.Y.Z can be listed with ``git log vX.Y.Z..`` for
    quick inspection.

2. Update the ``CHANGES.rst`` file to state the current date for this release
   and edit ``mutalyzer/__init__.py`` by updating `__date__` and removing the
   ``dev`` value from `__version_info__`.

   Commit and tag the version update::

       git commit -am 'Bump version to X.Y.Z'
       git tag -a 'vX.Y.Z'

3. Add a new entry at the top of the ``CHANGES.rst`` file like this::

       Version X.Y.Z+1
       ---------------

       Release date to be decided.

   Increment the patch version and add a ``dev`` value to `__version_info__`
   in ``mutalyzer/__init__.py`` and commit these changes::

       git commit -am 'Open development for X.Y.Z+1'

4. Push these commits to GitHub::

       git push github release-X.Y.Z -u

   And submit a pull request for this branch.

5. If everything looks ok and the pull request has been accepted you can push
   the tag and update the `release` branch. Beware to re-tag if the pull
   request was updated meanwhile. The working branch can be deleted.

   ::

       git push github +vX.Y.Z~0:refs/heads/release --tags
       git branch -d release-X.Y.Z

   That last push command might seem a bit cryptic (`it is explained here
   <http://stackoverflow.com/a/4061542>`_). It sets the remote branch
   `release` to whatever the tag `vX.Y.Z` points to and also pushes all tags.

   If the `release` branch is `protected
   <https://help.github.com/articles/about-protected-branches/>`_, updating it
   will be rejected until the required status checks pass (e.g., `Travis CI
   <http://travis-ci.org/>`_).


.. _SemVer: http://semver.org/
