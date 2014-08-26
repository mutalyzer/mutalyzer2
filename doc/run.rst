.. highlight:: none

.. _run:

Running Mutalyzer
=================

Please make sure Mutalyzer can find its configuration file, as detailed in
:ref:`config`.

Mutalyzer comes with a number of different interfaces, of which the website is
perhaps the main one. It can be started using a built-in test server that's
useful for development and debugging purposes like this::

    $ mutalyzer-website
     * Running on http://127.0.0.1:5000/

You can now point your webbrowser to the URL that is printed and see the
welcoming Mutalyzer homepage.

Likewise, the SOAP and HTTP/RPC+JSON webservices can be started with the
``mutalyzer-service-json`` and ``mutalyzer-service-soap`` commands,
respectively.

For processing batch jobs, the batch processor must be running. This process
can be started from the command line and will keep running until it is stopped
by pressing Ctrl+C::

    $ mutalyzer-batch-processor
    ^Cmutalyzer-batch-processor: Hitting Ctrl+C again will terminate any running job!
    mutalyzer-batch-processor: Graceful shutdown

The built-in test servers won't get you far in production, though, and there
are many other possibilities for deploying Mutalyzer using WSGI. This topic is
discussed in :ref:`deploy`.
