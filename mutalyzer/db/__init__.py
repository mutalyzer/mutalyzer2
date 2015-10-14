"""
This package contains database models and a scoped session definition, all
using SQLAlchemy.
"""


from __future__ import unicode_literals

import sqlalchemy
from sqlalchemy.engine.url import make_url
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy.pool import StaticPool

from mutalyzer.config import settings


class SessionFactory(sessionmaker):
    """
    Session factory that configures the engine lazily at first use with the
    current settings.DATABASE_URI.
    """
    def __call__(self, **local_kw):
        if self.kw['bind'] is None and 'bind' not in local_kw:
            self.kw['bind'] = create_engine()
        return super(SessionFactory, self).__call__(**local_kw)


def create_engine():
    """
    Create an SQLAlchemy connection engine from the current configuration.
    """
    if not settings.DATABASE_URI:
        # Just return silently when no database is configured (this function
        # may still be called via the configuration settings hook). Of course
        # actually using the database will fail.
        return

    url = make_url(settings.DATABASE_URI)
    options = {}

    if settings.DEBUG:
        options.update(echo=True)

    if url.drivername == 'sqlite' and url.database in (None, '', ':memory:'):
        # SQLite in-memory database are created per connection, so we need a
        # singleton pool if we want to see the same database across threads,
        # web requests, etcetera.
        options.update(
            connect_args={'check_same_thread': False},
            poolclass=StaticPool)

        engine = sqlalchemy.create_engine(url, **options)

        # For convenience, we also create tables if we're using an SQLite
        # in-memory database. By definition they won't yet exist.
        Base.metadata.create_all(engine)
        return engine

    return sqlalchemy.create_engine(url, **options)


def configure_session(uri):
    """
    (Re)configure the session by closing the existing session if it exists and
    loading the current configuration for use by future sessions.
    """
    global session_factory, session
    session.remove()
    session_factory.configure(bind=create_engine())


# Reconfigure the session if database configuration is updated.
settings.on_update(configure_session, 'DATABASE_URI')


# Sessions are automatically created where needed and are scoped by thread.
session_factory = SessionFactory()

#: Global scoped :class:`sqlalchemy.orm.session.Session` instance. Use this
#: for all database communication, except for querying models. For the latter,
#: each model has a `query` property that is a
#: :class:`sqlalchemy.orm.query.Query` object against the model and the
#: current `Session` when called.
session = scoped_session(session_factory)


#: Base class to use for our models.
Base = declarative_base()
Base.query = session.query_property()
