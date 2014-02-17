Changelog
=========

This is a record of changes made between each Mutalyzer release.


Version 2.0.beta-31
-------------------

Release date to be decided.


Version 2.0.beta-30
-------------------

Released on February 18th 2014.

- Handle NCBI Entrez response validation errors (fixes, among other things,
  `LOVD#29 <https://humgenprojects.lumc.nl/trac/LOVD3/ticket/29>`_).
- Loosen error severity when CDS cannot be translated.
- Mutalyzer development migrated from Subversion to Git for version control.


Version 2.0.beta-29
-------------------

Released on October 11th 2013.

- Add Jonathan Vis attribution and COMMIT logo to about page.


Version 2.0.beta-28
-------------------

Released on September 18th 2013.

- Enable the HTTP/RPC+JSON web service to be used with POST requests.


Version 2.0.beta-27
-------------------

Released on June 18th 2013.

- Fix caching transcript-protein links from NCBI, reducing impact of NCBI
  communication problems.


Version 2.0.beta-26
-------------------

Released on April 9th 2013.

- Added mm10 (Mouse) transcript mappings to position converter.
- LRG parser updated to LRG 1.7 schema (`#127
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/127>`_).


Version 2.0.beta-25
-------------------

Released on March 25th 2013.

- Detect incorrect exon annotation in transcript references.
- Move documentation to Trac.
- Exon table is included in `runMutalyzer` webservice results.
- Temporarily disable frameshift detection in experimental description
  extractor (`#124
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/124>`_).
- Allow selectors on transcript references in position converter.
- Syntax checker now supports protein level variant descriptions.


Version 2.0.beta-24
-------------------

Released on December 10th 2012.

- Rename some warning codes (webservice API) (`#98
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/98>`_).
- Variants on mtDNA in position converter.


Version 2.0.beta-23
-------------------

Released on November 8th 2012.

No user-visible changes.


Version 2.0.beta-22
-------------------

Released on November 2nd 2012.

- Submitting batch jobs via the web services (`#115
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/115>`_).
- Allow for leading whitespace in batch job input (`#107
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/107>`_).
- New `descriptionExtract` webservice function.
- Name checker now includes description extractor output as an experimental
  service.
- Slice chromosome by gene name in reference file loader is now case
  insensitive (`#118
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/118>`_).
- Warn on missing positioning scheme (`#114
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/114>`_).


Version 2.0.beta-21
-------------------

Released on July 23rd 2012.

- Support compound variants in position converter.
- Support non-coding transcripts in position converter (`#102
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/102>`_).
- Move to new RPC library version, causing slight change in HTTP/RPC+JSON
  webservice output (more wrappers around output), but fixes #104.
- Fix position converter for delins with explicit deleted sequence.
- Fix description update from Version 2.0.beta-20 to use- notation instead of
  counting.


Version 2.0.beta-20
-------------------

Released on July 21st 2012.

- Disabled the ``-u`` and ``+d`` convention in favour of the official HGVS
  recommendations.


Version 2.0.beta-19
-------------------

Released on June 21st 2012.

- Fix crash on inversions (`#99
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/99>`_).


Version 2.0.beta-18
-------------------

Released on June 7th 2012.

- Moved from soaplib to rpclib for webservices (`#66
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/66>`_).
- Added HTTP/RPC+JSON webservice (`#18
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/18>`_).
- Fixed name checker errors in some adjacent variants (`#83
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/83>`_).
- Name checker form now uses GET requests to support easier linking to result
  pages.
- You can now specify chromosomes by name in the reference file loader (`#92
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/92>`_).
- Made batch daemon not crash on MySQL restarts (`#91
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/91>`_).
- Position converter now detects incorrect order in position ranges (`#95
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/95>`_).
- Added NBIC logo to 'about' page.


Version 2.0.beta-17
-------------------

Released on April 2nd 2012.

- Fixed crossmapping bug for some transcripts.
- Fixes for NCBI Entrez EFetch Version 2.0 release.
- Better chromosomal variant descriptions.
- Various smaller features and bugfixes.


Version 2.0.beta-16
-------------------

Released on March 1st 2012.

- Fixed position converter mapping info for some transcripts.
- Fixed deletion with deleted sequence length as argument.


Version 2.0.beta-15
-------------------

Released on February 20th 2012.

- Added 'Description Extractor' (see the main menu).
- Fixes for NCBI Entrez EFetch Version 2.0 release.
- Added chromosomal positions to `getTranscriptsAndInfo` webservice.
- Fixed chromosome slicing on reverse complement
- Fixed describing NOP variants with ``=``.
- Added Reference sequence info in `runMutalyzer` SOAP function response.
- Fixed mapping info for genes mapped to more than one chromosome.
- Various smaller features and bugfixes.


Version 2.0.beta-14
-------------------

Released on January 26th 2012.

- Added a SOAP service `getTranscriptsMapping`.
- Various smaller features and bugfixes.


Version 2.0.beta-13
-------------------

Released on January 25th 2012.

- Accept EX positioning scheme.
- Fix handling of LRG reference sequences.
- Various smaller features and bugfixes.


Version 2.0.beta-12
-------------------

Released on November 25th 2011.

- Accept plasmid reference sequences.
- View variant position in UCSC Genome Browser (only for transcript
  references).
- Retry querying dbSNP if it does not respond the first time.
- Support reference GenBank files built from contigs.
- Add optional argument to SOAP service `numberConversion` to map chromosomal
  locations to any gene.
- Various smaller features and bugfixes.


Version 2.0.beta-11
-------------------

Released on September 30st 2011.

- Major code refactoring:

  - Mutalyzer is now structured as a proper Python package.
  - Reworked installation and upgrade procedure.
  - Remote installation using Fabric.
  - Batch scheduler is now a proper system daemon.
  - Use mod_wsgi (with web.py) instead of the deprecated mod_python.
  - Added a lot of internal documentation.
  - Introduce unit tests.
  - Handle deletions of entire exons.
  - Added a SOAP service `info`.
  - Handle unknown (fuzzy) intronic positions.
  - Automatic synchronization of database and cache between Mutalyzer
    installations.
  - Use NCBI instead of UCSC for transcript mapping info.
  - Added a SOAP service `getdbSNPDescriptions`.
  - Moved Trac and Subversion repository to new server.
  - Implement HTTP HEAD method for ``/Reference/*`` locations.

- Added a SOAP service `ping`.
- Added an optional versions parameter to the SOAP service `getTranscripts`.
- Various smaller features and bugfixes.


Version 2.0.beta-10
-------------------

Released on July 21st 2011.

- Greatly reduce runtime for large batch jobs.


Version 2.0.beta-9
------------------

Released on June 27th 2011.

- Reworked the calculation of new splice site positions.
- Optionally restrict SOAP service `getTranscriptsAndInfo` transcripts to a
  gene.
- Add raw variants to SOAP service `runMutalyzer` results.
- Provide webservice client examples.
- Various smaller features and bugfixes.


Older versions
--------------

The first lines of code for Mutalyzer 2.0 where written July 28th 2009, and
version 2.0.beta-8 was released on January 31st 2011. As far as Mutalyzer 1 is
concerned, archaeology is not really our field of research.
