.. highlight:: none

.. _download:

Downloading Mutalyzer
=====================

The Mutalyzer source code is `hosted on the LUMC GitLab server
<https://git.lumc.nl/mutalyzer/mutalyzer>`_. The recommended way to get the
Mutalyzer source code is by cloning the `Git`_ repository::

    git clone https://git.lumc.nl/mutalyzer/mutalyzer.git

This will give you the current development version. See below for working with
other versions.


Release versions
----------------

All Mutalyzer releases are tagged. You can run ``git tag`` to list the
available tags. Use ``git checkout <tag>`` to switch to a certain release. For
example, to swith to the `2.0.0` release::

    git checkout v2.0.0

Alternatively, you can checkout the `release` branch. This branch always
points to the latest Mutalyzer release.


Archive downloads
-----------------

If for whatever reason you don't want to use Git, you can download the source
code directly as a zip archive or tarball. The current development version can
be found from the `project homepage
<https://git.lumc.nl/mutalyzer/mutalyzer>`_. Archive downloads for release
versions can be found on the `tags page
<https://git.lumc.nl/mutalyzer/mutalyzer/tags>`_.


.. _Git: http://git-scm.com/