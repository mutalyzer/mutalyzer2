"""
Mutalyzer website interface using the Flask framework.
"""


import pkg_resources

from flask import Flask

from mutalyzer.config import settings
from mutalyzer.db import session


# Todo: Perhaps we also need this for the RPC services?
class ReverseProxied(object):
    """
    Wrap the application in this middleware and configure the front-end server
    to add these headers, to let you quietly bind this to a URL other than /
    and to an HTTP scheme that is different than what is used locally.

    Example for nginx::

        location /myprefix {
            proxy_pass http://192.168.0.1:5001;
            proxy_set_header Host $host;
            proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
            proxy_set_header X-Scheme $scheme;
            proxy_set_header X-Script-Name /myprefix;
        }

    `Flask Snippet <http://flask.pocoo.org/snippets/35/>`_ from Peter Hansen.
    """
    def __init__(self, app):
        self.app = app

    def __call__(self, environ, start_response):
        script_name = environ.get('HTTP_X_SCRIPT_NAME', '')
        if script_name:
            environ['SCRIPT_NAME'] = script_name
            path_info = environ['PATH_INFO']
            if path_info.startswith(script_name):
                environ['PATH_INFO'] = path_info[len(script_name):]

        scheme = environ.get('HTTP_X_SCHEME', '')
        if scheme:
            environ['wsgi.url_scheme'] = scheme
        return self.app(environ, start_response)


def create_app():
    """
    Create a Flask instance for Mutalyzer.
    """
    template_folder = pkg_resources.resource_filename(
        'mutalyzer', 'website/templates')
    static_folder = pkg_resources.resource_filename(
        'mutalyzer', 'website/templates/static')

    app = Flask('mutalyzer',
                template_folder=template_folder, static_folder=static_folder)

    app.config.update(DEBUG=settings.DEBUG,
                      TESTING=settings.TESTING,
                      MAX_CONTENT_LENGTH=settings.MAX_FILE_SIZE)
    from mutalyzer.website.views import website
    app.register_blueprint(website)

    @app.teardown_appcontext
    def shutdown_session(exception=None):
        session.remove()

    return app


def create_reverse_proxied_app():
    """
    Create a Flask instance for Mutalyzer running behind a reverse proxy.

    See :func:`create_app` and :class:`ReverseProxied`.
    """
    app = create_app()
    app.wsgi_app = ReverseProxied(app.wsgi_app)
    return app
