API reference
=============

.. note:: This reference is incomplete. In particular, modules without
          reStructuredText formatted docstrings are omitted.


.. mutalyzer
   ---------

   .. automodule:: mutalyzer
      :members: NOMENCLATURE_VERSION_INFO, NOMENCLATURE_VERSION, COPYRIGHT_YEARS, SOAP_NAMESPACE


mutalyzer.announce
------------------

.. automodule:: mutalyzer.announce
   :members:


mutalyzer.config
----------------

.. automodule:: mutalyzer.config
   :members:


.. mutalyzer.config.default_settings
   ---------------------------------

   .. automodule:: mutalyzer.config.default_settings
      :members:


.. mutalyzer.Crossmap
   ------------------

   .. automodule:: mutalyzer.Crossmap
      :members:


mutalyzer.db
------------

.. automodule:: mutalyzer.db
   :members:


mutalyzer.db.models
-------------------

.. I have no idea why, but somehow :members: includes `or_` and `relationship`
   from SQLAlchemy (which are imported in the module).

.. automodule:: mutalyzer.db.models
   :members:
   :exclude-members: or_, relationship


mutalyzer.db.queries
--------------------

.. automodule:: mutalyzer.db.queries
   :members:
   :exclude-members: and_, or_


.. mutalyzer.describe
   ------------------

   .. automodule:: mutalyzer.describe
      :members:


.. mutalyzer.entrypoints
   ---------------------

   .. automodule:: mutalyzer.entrypoints
      :members:


mutalyzer.entrypoints.admin
---------------------------

.. automodule:: mutalyzer.entrypoints.admin
   :members:


mutalyzer.entrypoints.batch_processor
-------------------------------------

.. automodule:: mutalyzer.entrypoints.batch_processor
   :members:


mutalyzer.entrypoints.mutalyzer
-------------------------------

.. automodule:: mutalyzer.entrypoints.mutalyzer
   :members:


mutalyzer.entrypoints.service_json
----------------------------------

.. automodule:: mutalyzer.entrypoints.service_json
   :members:


mutalyzer.entrypoints.service_soap
----------------------------------

.. automodule:: mutalyzer.entrypoints.service_soap
   :members:


mutalyzer.entrypoints.website
-----------------------------

.. automodule:: mutalyzer.entrypoints.website
   :members:


.. mutalyzer.File
   --------------

   .. automodule:: mutalyzer.File
      :members:


.. mutalyzer.GenRecord
   -------------------

   .. automodule:: mutalyzer.GenRecord
      :members:


.. mutalyzer.grammar
   -----------------

   .. automodule:: mutalyzer.grammar
      :members:


.. mutalyzer.mapping
   -----------------

   .. automodule:: mutalyzer.mapping
      :members:


.. mutalyzer.models
   ----------------

   .. automodule:: mutalyzer.models
      :members:


.. mutalyzer.mutator
   -----------------

   .. automodule:: mutalyzer.mutator
      :members:


.. mutalyzer.output
   ----------------

   .. automodule:: mutalyzer.output
      :members:


.. mutalyzer.parsers
   -----------------

   .. automodule:: mutalyzer.parsers
      :members:


.. mutalyzer.parsers.genbank
   -------------------------

   .. automodule:: mutalyzer.parsers.genbank
      :members:


.. mutalyzer.parsers.lrg
   ---------------------

   .. automodule:: mutalyzer.parsers.lrg
      :members:


mutalyzer.redisclient
---------------------

.. automodule:: mutalyzer.redisclient
   :members:


.. mutalyzer.Retriever
   -------------------

   .. automodule:: mutalyzer.Retriever
      :members:


.. mutalyzer.Scheduler
   -------------------

   .. automodule:: mutalyzer.Scheduler
      :members:


.. mutalyzer.services
   ------------------

   .. automodule:: mutalyzer.services
      :members:


mutalyzer.services.json
-----------------------

.. automodule:: mutalyzer.services.json
   :members:


.. mutalyzer.services.rpc
   ----------------------

   .. automodule:: mutalyzer.services.rpc
      :members:


mutalyzer.services.soap
-----------------------

.. automodule:: mutalyzer.services.soap
   :members:


mutalyzer.stats
---------------

.. automodule:: mutalyzer.stats
   :members:


mutalyzer.sync
--------------

.. automodule:: mutalyzer.sync
   :members:


.. mutalyzer.util
   --------------

   .. automodule:: mutalyzer.util
      :members:


.. mutalyzer.variantchecker
   ------------------------

   .. automodule:: mutalyzer.variantchecker
      :members:


mutalyzer.website
-----------------

.. automodule:: mutalyzer.website
   :members:


mutalyzer.website.views
-----------------------

.. automodule:: mutalyzer.website.views
   :members:
