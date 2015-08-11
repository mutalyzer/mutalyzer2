.. highlight:: none

.. _deploy:

Deploying Mutalyzer in production
=================================

The previous sections discussed managing a Mutalyzer installation with a focus
on a development environment. There are a number of additional things you will
want to consider when deploying Mutalyzer in a production environment, mainly
concerning security and performance.

Usually you'll at least want to use a well-performing WSGI application server
for the website and SOAP and HTTP/RPC+JSON webservices. There are many options
here, ranging from Apache's `mod_wsgi`_ to `uWSGI`_ to standalone WSGI
containers such as `Gunicorn`_.

Below we briefly describe our recommended setup for a production environment
using Gunicorn, nginx and Supervisor.


Configuration settings
----------------------

Todo: Link to the description of these configuration settings.

It is recommended to at least set the following configuration settings:

- DEBUG
- EMAIL
- CACHE_DIR
- SOAP_WSDL_URL
- JSON_ROOT_URL


WSGI application server: Gunicorn
---------------------------------

`Gunicorn`_ is a well-perfoming Python WSGI HTTP Server. Being a Python
application, it can be installed in the Mutalyzer virtual environment with
``pip install gunicorn``.

Many configuration settings are available for Gunicorn and we recommend to use
a configuration file per WSGI application. For example, the following
configuration can be stored in ``website.conf``:

.. code-block:: ini

    workers = 4
    max_requests = 1000
    timeout = 600
    bind = 'unix:/opt/mutalyzer/run/website.sock'

This will bind the Gunicorn server to a unix socket (which we can later use
from nginx) and run with 4 worker processes. To serve the Mutalyzer website
with this configuration, run the following::

    $ gunicorn -c website.conf mutalyzer.entrypoints.website

This uses the WSGI application object exported by the
`mutalyzer.entrypoints.website` module. Likewise, the SOAP and HTTP/RPC+JSON
webservices have WSGI application objects exported by the
`mutalyzer.entrypoints.service_soap` and `mutalyzer.entrypoints.service_json`
modules.


Web server: nginx
-----------------

It is usually a good idea to use a separate webserver in front of the WSGI
application servers. We use `nginx`_ for this purpose and configure it to
server static files directly and act as a reverse proxy for the WSGI
applications.

For example, to serve the website from the root path and the HTTP/RPC+JSON
webservice from the ``/json`` path, an nginx configuration similar to the
following can be used:

.. code-block:: nginx

    server {
        listen                     80;
        server_name                _;

        client_max_body_size       2G;
        keepalive_timeout          5;

        location /static/ {
            alias /opt/mutalyzer/static/;
            expires 30d;
            add_header Pragma public;
            add_header Cache-Control "public";
        }

        location / {
            root                   /usr/share/nginx/html;
            proxy_set_header       X-Forwarded-For $proxy_add_x_forwarded_for;
            proxy_set_header       X-Real-IP $remote_addr;
            proxy_set_header       X-Scheme $scheme;
            proxy_set_header       Host $http_host;
            proxy_redirect         off;
            proxy_read_timeout     600;
            proxy_pass             http://website;
        }

        location /json {
            root                   /usr/share/nginx/html;
            proxy_set_header       X-Forwarded-For $proxy_add_x_forwarded_for;
            proxy_set_header       X-Real-IP $remote_addr;
            proxy_set_header       X-Scheme $scheme;
            proxy_set_header       Host $http_host;
            proxy_redirect         off;
            proxy_read_timeout     600;
            proxy_pass             http://service-json;
        }
    }

    upstream website {
        server                     unix:/opt/mutalyzer/run/website.sock fail_timeout=0;
    }

    upstream service-json {
        server                     unix:/opt/mutalyzer/run/service-json.sock fail_timeout=0;
    }


Process control: Supervisor
---------------------------

For managing the different WSGI application servers and Mutalyzer batch
processor, Supervisor can be used. Supervisor is usually started from the init
system and controls programs and program groups. For example, it can
automatically restart a program if it crashed for some reason.

The following is an example Supervisor configuration defining a Mutalyzer
group consisting of the batch processor and a Gunicorn process for the
website:

.. code-block:: ini

    [group:mutalyzer]
    programs=batch-processor,website

    [program:batch-processor]
    command=mutalyzer-batch-processor
    autorestart=true
    environment=MUTALYZER_SETTINGS="/opt/mutalyzer/conf/settings.py"

    [program:website]
    command=gunicorn -c /opt/mutalyzer/conf/website.conf mutalyzer.entrypoints.website
    autorestart=true
    environment=MUTALYZER_SETTINGS="/opt/mutalyzer/conf/settings.py"


Automated deployment with Ansible
---------------------------------

Deployments of complete production environments are often complex and
repetitive. Therefore, manual deployments are inefficient and
error-prone. Several systems exist to automate this, such as `Puppet`_,
`Chef`_, and `Ansible`_.

An automated `deployment of Mutalyzer with Ansible
<https://github.com/mutalyzer/ansible-role-mutalyzer>`_ is available on
GitHub. This includes installation of the website, SOAP and HTTP/RPC+JSON
webservices, and the batch processor, similar to the setup described above.


.. _Ansible: http://www.ansible.com/
.. _Chef: http://www.getchef.com/
.. _Gunicorn: http://gunicorn.org/
.. _mod_wsgi: https://code.google.com/p/modwsgi/
.. _nginx: http://nginx.org/
.. _Puppet: http://puppetlabs.com/
.. _uWSGI: http://uwsgi-docs.readthedocs.org/
