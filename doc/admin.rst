.. highlight:: none

.. _admin:

Administration
==============

For administrative tasks, Mutalyzer comes with the ``mutalyzer-admin`` command
line utility. Use its ``-h`` argument in combination with any subcommand for
detailed usage information, for example::

    $ mutalyzer-admin setup-database -h
    usage: mutalyzer-admin setup-database [-h] [--destructive] [-c ALEMBIC_CONFIG]

    Setup database tables (if they do not yet exist).

    optional arguments:
      -h, --help            show this help message and exit
      --destructive         delete any existing tables and data
      -c ALEMBIC_CONFIG, --alembic-config ALEMBIC_CONFIG
                            path to Alembic configuration file

    If Alembic config is given (--alembic-config), this also prepares the database
    for future migrations with Alembic (recommended).


Managing genome assemblies
--------------------------

Mutalyzer can be loaded with any number of genome assemblies. Each assembly
includes a list of chromosomes. To list the currently loaded genome
assemblies::

    $ mutalyzer-admin assemblies list
    GRCh37 (hg19), Homo sapiens (9606)
    GRCm38 (mm10), Mus musculus (10090)

Loading a new genome assembly is done with information in a JSON file, of
which there are some examples in the Mutalyzer source tree under the
``extras/assemblies`` directory, for example::

    $ mutalyzer-admin assemblies add extras/assemblies/GRCh37.json

For any genome assembly, transcript mappings can be imported. These include
genomic coordinate mappings of their CDS and exons. Currently, three sources
of transcript mappings are supported.

.. note:: The following ``mutalyzer-admin assemblies`` subcommands all accept
          an optional ``--assembly`` argument to specify the genome assembly
          to import to.


Import mappings from an NCBI mapview file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The NCBI provides FTP downloads of transcript mappings for a large number of
genome assemblies as used by their `Map Viewer
<http://www.ncbi.nlm.nih.gov/mapview/>`_ service. These can be imported with
``mutalyzer-admin``, but only after sorting by the *feature_id* and
*chromosome* columns.

For example, to import transcript mappings for the GRCh37 assembly, run the
following::

    $ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/MapView/Homo_sapiens/sequence/ANNOTATION_RELEASE.105/initial_release/seq_gene.md.gz
    $ zcat seq_gene.md.gz | sort -t $'\t' -k 11,11 -k 2,2 > seq_gene.sorted.md
    $ mutalyzer-admin assemblies import-mapview seq_gene.sorted.md 'GRCh37.p13-Primary Assembly'

.. note:: The last argument, ``GRCh37.p13-Primary Assembly``, defines the group
          label to filter the file on. You would usually want to include it.

Examples for other assemblies can be found in `this Gist
<https://gist.github.com/martijnvermaat/ce84945d05b4e42d3584>`_.


Import mappings from the UCSC Genome Browser MySQL database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Transcript mappings from the `UCSC Genome Browser MySQL database
<https://genome.ucsc.edu/goldenPath/help/mysql.html>`_ can be imported on a
per-gene basis. This is useful when the NCBI mappings do not (yet) include a
certain gene or transcript.

For example, to import all TTN transcript mappings::

    $ mutalyzer-admin assemblies import-gene TTN

.. note:: This subcommand chooses the UCSC genome assembly by using the alias
          of the specified Mutalyzer genome assembly (`hg19` by default).


Import mappings from a reference file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For transcript mappings that are not available from our usual sources,
importing from a genomic reference is supported::

    $ mutalyzer-admin assemblies import-reference NC_012920.1

.. note:: Currently this subcommand is restricted to importing mtDNA
          transcripts, since it has the chromosome hard coded and only
          supports one exon per transcript.


Showing announcements to users
------------------------------

It is possible to define an announcement to be shown on the website
interface. For example, to display *Hello World!* with a link to the `GNU
Hello World! <http://www.gnu.org/fun/jokes/helloworld.html>`_ page::

    $ mutalyzer-admin announcement set 'Hello World!' \
        --url http://www.gnu.org/fun/jokes/helloworld.html

To remove the announcement, use ``unset``::

    $ mutalyzer-admin announcement unset


Synchronizing the cache with other installations
------------------------------------------------

Using the ``sync-cache`` subcommand, the reference file cache of a remote
Mutalyzer installation can be queried for new entries which are then retrieved
and added to the local cache.

The primary purpose for this is synchronizing reference files loaded by users
with the reference file loader between different servers. These reference
files are assigned a unique accession number (starting with ``UD_``) upon
creation, which is at that point unknown to any other Mutalyzer server.

For example, to synchronize the local reference file cache with the `primary
Mutalyzer server <https://mutalyzer.nl/>`_::

    $ mutalyzer-admin sync-cache 'https://mutalyzer.nl/services/?wsdl' \
        'https://mutalyzer.nl/Reference/{file}'


Mutalyzer database setup
------------------------

After installation, a database needs to be setup for Mutalyzer to run (see
:ref:`install-setup`)::

    $ mutalyzer-admin setup-database --alembic-config migrations/alembic.ini

The ``--alembic-config`` argument points to the ``alembic.ini`` file in the
Mutalyzer source tree and it enables initialization of database migration
management. It is recommended to include it, but you don't need it if you
don't plan to ever upgrade your Mutalyzer installation.

This subcommand also takes an optional ``--destructive`` argument, which can
be used to remove any existing database content.
