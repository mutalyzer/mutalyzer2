"""
Mutalyzer website interface using the Flask framework.
"""


from __future__ import unicode_literals

import logging
import os
import pkg_resources
import urlparse

from flask import Flask

from mutalyzer.config import settings
from mutalyzer.db import session


def create_app():
    """
    Create a Flask instance for Mutalyzer.
    """
    template_folder = pkg_resources.resource_filename(
        'mutalyzer', 'website/templates')
    static_folder = pkg_resources.resource_filename(
        'mutalyzer', 'website/templates/static')

    app = Flask('mutalyzer',
                template_folder=os.path.abspath(template_folder),
                static_folder=os.path.abspath(static_folder))

    app.config.update(DEBUG=settings.DEBUG,
                      TESTING=settings.TESTING,
                      MAX_CONTENT_LENGTH=settings.MAX_FILE_SIZE)
    from mutalyzer.website.views import website
    app.register_blueprint(website)

    @app.before_first_request
    def setup_logging():
        if not settings.DEBUG:
            # In production mode, log errors to standard error.
            app.logger.addHandler(logging.StreamHandler())
            app.logger.setLevel(logging.ERROR)

    @app.teardown_appcontext
    def shutdown_session(exception=None):
        session.remove()

    return app


def url_for(endpoint, **values):
    """
    Generates a URL to the given website endpoint.

    Like :func:`Flask.url_for`, but for when you don't have an application or
    request context.

    Note that the generated URL will be based on the `WEBSITE_ROOT_URL`
    configuration setting or `http://localhost` if not set.

    :arg str endpoint: The endpoint of the URL (name of the function).
    :arg str values: The variable arguments of the URL rule.
    """
    root = urlparse.urlsplit(settings.WEBSITE_ROOT_URL or 'http://localhost')
    url_map = create_app().url_map.bind(root.netloc, root.path or '/',
                                        url_scheme=root.scheme)
    return url_map.build('website.%s' % endpoint, values, force_external=True)
