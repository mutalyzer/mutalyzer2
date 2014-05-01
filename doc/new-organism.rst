Adding a new organism to Mutalyzer
==================================


Introduction
------------

In this document, we describe what is needed for Mutalyzer to support new
organisms. For each functionality, we give a list of requirements and an
estimate of the amount of time needed for implementation.


Position converter
------------------

To use the position converter, we need a mapping database for the genomic
reference sequence that is used. The database should contain the following
information per gene:

- Genomic location (chromosome, transcription start, transcription end).
- Gene model (CDS start, CDS stop, all splice sites).

The transcription start and end may be missing, but this is not recommended.


Implementation time
^^^^^^^^^^^^^^^^^^^

Depending on the format of the database, making this functionality available
is relatively straightforward. If the database is stored in a structured
format (CSV or something that can be converted to CSV automatically),
importing should take no more than two working days.


Name checker
------------

To use the name checker, a fully annotated GenBank record should be available
for every genomic location.


Public genome build
^^^^^^^^^^^^^^^^^^^

If the genome build of the organism in question is public, as is the case for
some model organisms, the GenBank records can be retrieved from the NCBI. In
this case, no special effort is required.

If the data is public, but is not yet available as a genome build at the NCBI,
the mapping database can be converted to the format required by the NCBI and
uploaded. The conversion should take no more than two working days, but the
time it will take before the build is available depends on the NCBI.

Non public genome build
^^^^^^^^^^^^^^^^^^^^^^^

If it is not possible to submit the generated GenBank reference files to the
NCBI, for commercial reasons for example, there is the option to offer a
stand-alone version of Mutalyzer containing the reference files. This version
should be hosted at the client side. Preparing such an installation will
require no more than three working days.


Additional notes
----------------

There is currently no program available for the generation of the GenBank
reference files. Although this is a one time effort, we estimate that the
development of such a program requires three weeks.

There is also no program available for the conversion of a mapping database to
the format the NCBI expects. The specifications are also unknown. Once we have
the specifications, we estimate that the development of a conversion tool will
take three weeks.

A stand-alone version of Mutalyzer will require some expertise to set up and
maintain at the client side. Although we currently provide no support or
maintenance, we are thinking about a model for this.
