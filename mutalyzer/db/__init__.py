"""
This package contains database models and a scoped session definition, all
using SQLAlchemy.
"""


from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import scoped_session, sessionmaker


engine = create_engine('sqlite:////tmp/test.db') #, echo=True)

session_factory = sessionmaker(bind=engine)

session = scoped_session(session_factory)

Base = declarative_base()
Base.query = session.query_property()


def create_database():
    Base.metadata.drop_all(engine)
    Base.metadata.create_all(engine)

    # if using alembic:
    #from alembic.config import Config
    #from alembic import command
    #alembic_cfg = Config("alembic.ini")
    #command.stamp(alembic_cfg, "head")
