"""
Mutalyzer website interface using the Flask framework.
"""


import logging
import os
import pkg_resources

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
