Changelog
=========

This is a record of changes made between each Mutalyzer release.


Version 2.0.14
--------------

Release date to be decided.

- E-mail address is now optional in batch job website interface (`#119
  <https://github.com/mutalyzer/mutalyzer/pull/119>`_).
- Use Mailcheck in the batch jobs form (`#107
  <https://github.com/mutalyzer/mutalyzer/pull/107>`_).
- Optional email for batch jobs from webservices (`#104
  <https://github.com/mutalyzer/mutalyzer/pull/104>`_).
- Process batch jobs grouped by email address (`#101
  <https://github.com/mutalyzer/mutalyzer/pull/101>`_).
- Use interval binning scheme on transcript mappings (`#100
  <https://github.com/mutalyzer/mutalyzer/pull/100>`_).
- Back translator interface (`#74
  <https://github.com/mutalyzer/mutalyzer/pull/74>`_).
- Transcript-protein links are now cached in Redis (`#94
  <https://github.com/mutalyzer/mutalyzer/pull/94>`_).
- Added `Ter` as a valid amino acid in the HGVS grammar (`#90
  <https://github.com/mutalyzer/mutalyzer/pull/90>`_).
- Refactoring of unit tests (`#88
  <https://github.com/mutalyzer/mutalyzer/pull/88>`_, `#89
  <https://github.com/mutalyzer/mutalyzer/pull/89>`_).


Version 2.0.13
--------------

Released on October 1st 2015.

- Fix query bug in update transcript-protein links (`#85
  <https://github.com/mutalyzer/mutalyzer/pull/85>`_).


Version 2.0.12
--------------

Released on September 30th 2015.

- Fix off-by-one in slicing chromosome by gene name (`#79
  <https://github.com/mutalyzer/mutalyzer/pull/79>`_).
- Document scheme used for all positions and ranges (`#79
  <https://github.com/mutalyzer/mutalyzer/pull/79>`_).
- Show diff for variant protein from non-reference start codon (`#78
  <https://github.com/mutalyzer/mutalyzer/pull/78>`_).
- Visualise protein change, also with alternative start (`#72
  <https://github.com/mutalyzer/mutalyzer/pull/72>`_).
- Translate alternative start to M, also in variant (`#72
  <https://github.com/mutalyzer/mutalyzer/pull/72>`_).
- Added Baker's yeast (SacCer_Apr2011/sacCer3) assembly (`#73
  <https://github.com/mutalyzer/mutalyzer/pull/73>`_).


Version 2.0.11
--------------

Released on August 6th 2015.

- Fix bug in recognizing ``p.(=)`` (was reported as ``p.?``) (`#65
  <https://github.com/mutalyzer/mutalyzer/pull/65>`_).


Version 2.0.10
--------------

Released on July 21st 2015.

- Don't crash the position converter on transcript mappings containing no
  exons (`#63
  <https://github.com/mutalyzer/mutalyzer/pull/63>`_).
- Use the notation for an uncertain stop codon, e.g., ``p.(Gln730Profs*?)``
  instead of ``p.(Gln730Profs*96)`` when a variant results in a frame shift or
  extension and we don't see a new stop codon in the RNA (`#57
  <https://github.com/mutalyzer/mutalyzer/pull/57>`_).
- Added Dog (Broad CanFam3.1/canFam3) assembly for position converter (`#56
  <https://github.com/mutalyzer/mutalyzer/pull/56>`_).


Version 2.0.9
-------------

Released on July 9th 2015.

- Improvements in encoding detection of uploaded batch files (`#52
  <https://github.com/LUMC/mutalyzer/pull/52>`_, `#53
  <https://github.com/LUMC/mutalyzer/pull/53>`_).
- Usability improvements in reading DNA for description extractor (`#54
  <https://github.com/mutalyzer/mutalyzer/pull/54>`_).


Version 2.0.8
-------------

Released on May 27th 2015.

- The `HGVS variant description extractor
  <https://github.com/mutalyzer/description-extractor>`_ package is available
  through the Description Extractor interface (`#1
  <https://github.com/LUMC/mutalyzer/pull/1>`_).
- GitHub project moved from `LUMC/mutalyzer
  <https://github.com/LUMC/mutalyzer>`_ to `mutalyzer/mutalyzer
  <https://github.com/LUMC/mutalyzer>`_.


Version 2.0.7
-------------

Released on March 27th 2015.

- Return fault code and description on RPC service errors (`#31
  <https://github.com/LUMC/mutalyzer/pull/31>`_).
- Use esummary 2.0 response format (`#32
  <https://github.com/LUMC/mutalyzer/pull/32>`_).


Version 2.0.6
-------------

Released on February 10th 2015.

- Added `getGeneLocation` webservice method. Given a gene symbol and optional
  genome build, it returns the location of the gene (`#28
  <https://github.com/LUMC/mutalyzer/pull/28>`_).
- Discard incomplete genes in genbank reference files (`#26
  <https://github.com/LUMC/mutalyzer/pull/26>`_).


Version 2.0.5
-------------

Released on December 16th 2014.

- New website layout by `Landscape <http://wearelandscape.nl/>`_ (`GitLab!26
  <https://git.lumc.nl/mutalyzer/mutalyzer/merge_requests/26>`_).
- Source code moved from `GitLab <https://git.lumc.nl/mutalyzer/mutalyzer>`_
  to `GitHub <https://github.com/LUMC/mutalyzer>`_.
- Automated unit tests `on Travis CI <https://travis-ci.org/LUMC/mutalyzer>`_
  (`#16 <https://github.com/LUMC/mutalyzer/pull/16>`_).
- Developer documentation `hosted on Read the Docs
  <http://mutalyzer.readthedocs.org/en/latest/>`_ (`#17
  <https://github.com/LUMC/mutalyzer/pull/18>`_).


Version 2.0.4
-------------

Released on November 19th 2014.

- Many string encoding related fixes. Summarizing, Mutalyzer should now be
  completely aware of input and output string encodings, in all
  interfaces. Internally, all strings are unicode strings. This comes with
  minor changes in the webservice definitions, which most clients will
  probably not notice (`GitLab!25
  <https://git.lumc.nl/mutalyzer/mutalyzer/merge_requests/25>`_).
- Don't crash on mail errors in the batch scheduler (`GitLab#30
  <https://git.lumc.nl/mutalyzer/mutalyzer/issues/30>`_).
- Fix importing transcript mappings from UCSC database (`GitLab#9
  <https://git.lumc.nl/mutalyzer/mutalyzer/issues/9>`_).
- Rename GRCh36 to NCBI36 (`GitLab#8
  <https://git.lumc.nl/mutalyzer/mutalyzer/issues/8>`_).
- Updated all Python dependencies to their latest versions.
- Mutalyzer is now Open Source! Source code is available under the AGPL and
  documentation under the CC-by-sa license.


Version 2.0.3
-------------

Released on September 20th 2014.

- Fix several error cases in LOVD2 view on the name checker.


Version 2.0.2
-------------

Released on October 9th 2014.

- Fix incorrect GRCm38 chromosome accession number versions.
- Fix crash in position converter batch jobs.
- Upgrade the webservice library we use (Spyne, from 2.10.10 to 2.11.0). This
  potentially affects behaviour of both our SOAP and HTTP/RPC+JSON
  webservices, although our tests did not show any problems.


Version 2.0.1
-------------

Released on September 27th 2014.

- Fix POST requests to the HTTP/RPC+JSON webservice. This was a regression
  from version 2.0.beta-33. Thanks to Ken Doig for reporting the issue.


Version 2.0.0
-------------

Released on September 26th 2014.

This release does not bring many new features, but comes with significant
changes to the technical infrastructure. `GitLab!6
<https://git.lumc.nl/mutalyzer/mutalyzer/merge_requests/6>`_ tracks most of
this.

Some highlights especially users of the webservices should be aware of:

- HTTP/RPC+JSON webservice has changed response format (wrapper object
  removed). See below for an :ref:`example <changelog_200_example>`.
- No more plain HTTP access, only redirects to HTTPS.
- Many website entrypoints have changed URLs and form parameter names (the old
  ones have HTTP redirects).
- Removed old redirects from paths starting with ``/2.0/``.
- In maintenance mode, all requests get a *Service Temporarily Unavailable*
  response with status code 503.

Other changes:

- Upload a genbank file using the SOAP webservice (`uploadGenBankLocalFile`).
- Do not cleanup the cache during request handling (`GitLab#18
  <https://git.lumc.nl/mutalyzer/mutalyzer/issues/18>`_).
- Add GRCh38 (hg38) assembly (`GitLab!20
  <https://git.lumc.nl/mutalyzer/mutalyzer/merge_requests/20>`_).
- Move from nose to `pytest <http://pytest.org/>`_ for unit tests (`GitLab!23
  <https://git.lumc.nl/mutalyzer/mutalyzer/merge_requests/23>`_).
- Fix running Mutalyzer in a `virtual environment
  <http://virtualenv.readthedocs.org/>`_ and have an up-to-date
  ``requirements.txt`` for `pip <http://pip.readthedocs.org/>`_ (`GitLab!4
  <https://git.lumc.nl/mutalyzer/mutalyzer/merge_requests/4>`_).
- Switch from TAL to Jinja2 (`GitLab!3
  <https://git.lumc.nl/mutalyzer/mutalyzer/merge_requests/3>`_).
- Refactor user interfaces (`GitLab!5
  <https://git.lumc.nl/mutalyzer/mutalyzer/merge_requests/5>`_).
- Move from configobj to Python module based config (`GitLab!7
  <https://git.lumc.nl/mutalyzer/mutalyzer/merge_requests/7>`_).
- Use SQLAlchemy as ORM (`GitLab!8
  <https://git.lumc.nl/mutalyzer/mutalyzer/merge_requests/8>`_).
- Use Redis for stat counters (`GitLab!10
  <https://git.lumc.nl/mutalyzer/mutalyzer/merge_requests/10>`_).
- Port website from web.py to Flask (`GitLab!11
  <https://git.lumc.nl/mutalyzer/mutalyzer/merge_requests/11>`_).
- Isolated unit tests using fixtures and an in-memory database (`GitLab!12
  <https://git.lumc.nl/mutalyzer/mutalyzer/merge_requests/12>`_).
- Display announcement on website (`GitLab!14
  <https://git.lumc.nl/mutalyzer/mutalyzer/merge_requests/14>`_).
- Database migrations with Alembic (`GitLab!15
  <https://git.lumc.nl/mutalyzer/mutalyzer/merge_requests/15>`_).
- Update documentation and use Sphinx (`GitLab!16
  <https://git.lumc.nl/mutalyzer/mutalyzer/merge_requests/16>`_).
- Move to `semantic versioning <http://semver.org/>`_, starting with version
  2.0.0 (`GitLab!22
  <https://git.lumc.nl/mutalyzer/mutalyzer/merge_requests/22>`_).
- Add 404 not found page.
- Don't auto remove comma characters in syntax checker.
- Add a dash (``-``) as an allowed character in the gene name.
- Range, reverse complement range, and compound
  insertions/insertion-deletions.

.. _changelog_200_example:

The wrapper object has been removed from the HTTP/RPC+JSON webservice response
format. As an example, consider an old response format for the `checkSyntax`
method:

.. code-block:: json

    {
      "checkSyntaxResponse": {
        "checkSyntaxResult": {
          "valid": true,
          "messages": {
            "SoapMessage": []
          }
        }
      }
    }

The new response format is:

.. code-block:: json

    {
      "valid": true,
      "messages": []
    }


Version 2.0.beta-33
-------------------

Released on August 19th 2014.

- Link to `Upcoming server update
  <https://humgenprojects.lumc.nl/trac/mutalyzer/wiki/News/2014-08-19-upcoming-server-update>`_
  announcement.


Version 2.0.beta-32
-------------------

Released on June 26th 2014.

- Link to `Visual interface for Variant Description Extractor
  <https://humgenprojects.lumc.nl/trac/mutalyzer/wiki/News/2014-06-26-visual-interface>`_
  announcement.


Version 2.0.beta-31
-------------------

Released on March 27th 2014.

- Due to incorrect interpretation, temporarily only support one CDS per
  transcript (ignore all others) in LRG.
- Due to incorrect interpretation, temporarily ignore transcripts without a
  fixed id.


Version 2.0.beta-30
-------------------

Released on February 18th 2014.

- Handle NCBI Entrez response validation errors (fixes, among other things,
  `LOVD Trac#29 <https://humgenprojects.lumc.nl/trac/LOVD3/ticket/29>`_).
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
- LRG parser updated to LRG 1.7 schema (`Trac#127
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/127>`_).


Version 2.0.beta-25
-------------------

Released on March 25th 2013.

- Detect incorrect exon annotation in transcript references.
- Move documentation to Trac.
- Exon table is included in `runMutalyzer` webservice results.
- Temporarily disable frameshift detection in experimental description
  extractor (`Trac#124
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/124>`_).
- Allow selectors on transcript references in position converter.
- Syntax checker now supports protein level variant descriptions.


Version 2.0.beta-24
-------------------

Released on December 10th 2012.

- Rename some warning codes (webservice API) (`Trac#98
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/98>`_).
- Variants on mtDNA in position converter.


Version 2.0.beta-23
-------------------

Released on November 8th 2012.

No user-visible changes.


Version 2.0.beta-22
-------------------

Released on November 2nd 2012.

- Submitting batch jobs via the web services (`Trac#115
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/115>`_).
- Allow for leading whitespace in batch job input (`Trac#107
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/107>`_).
- New `descriptionExtract` webservice function.
- Name checker now includes description extractor output as an experimental
  service.
- Slice chromosome by gene name in reference file loader is now case
  insensitive (`Trac#118
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/118>`_).
- Warn on missing positioning scheme (`Trac#114
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/114>`_).


Version 2.0.beta-21
-------------------

Released on July 23rd 2012.

- Support compound variants in position converter.
- Support non-coding transcripts in position converter (`Trac#102
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/102>`_).
- Move to new RPC library version, causing slight change in HTTP/RPC+JSON
  webservice output (more wrappers around output), but fixes `Trac#104
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/104>`_.
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

- Fix crash on inversions (`Trac#99
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/99>`_).


Version 2.0.beta-18
-------------------

Released on June 7th 2012.

- Moved from soaplib to rpclib for webservices (`Trac#66
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/66>`_).
- Added HTTP/RPC+JSON webservice (`Trac#18
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/18>`_).
- Fixed name checker errors in some adjacent variants (`Trac#83
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/83>`_).
- Name checker form now uses GET requests to support easier linking to result
  pages.
- You can now specify chromosomes by name in the reference file loader
  (`Trac#92 <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/92>`_).
- Made batch daemon not crash on MySQL restarts (`Trac#91
  <https://humgenprojects.lumc.nl/trac/mutalyzer/ticket/91>`_).
- Position converter now detects incorrect order in position ranges (`Trac#95
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

The first lines of code for Mutalyzer 2.0 were written July 28th 2009, and
version 2.0.beta-8 was released on January 31st 2011. As far as Mutalyzer 1 is
concerned, archaeology is not really our field of research.
