"""
Fabric fabfile for Mutalyzer.

Notice: The definitions in this file are quite specific to the standard
Mutalyzer environment. This consists of a Debian stable (Squeeze) system with
Apache and Mutalyzer using its mod_wsgi module. Debian conventions are used
throughout. See the README file for more information.

To do a deployment on a server with an existing configured Mutalyzer
installation:

  $ fab deploy -H server1.mutalyzer.nl

For a fresh deployment on a new server:

  $ fab deploy:boostrap=yes -H server1.mutalyzer.nl
"""


from fabric.api import *


def deploy(bootstrap='no'):
    """
    Deploy Mutalyzer on the remote host.

    Create a source distribution, transfer it to the remote host, and install
    from there. After installation, we restart Apache and the Mutalyzer batch
    daemon.

    Additionally, if bootstrap=yes, install all dependencies before Mutalyzer
    installation, and bootstrap the Mutalyzer configuration afterwards (i.e.
    create and fill database, add cron script, create cache directory, etc).
    """
    # Currently, Fabric only supports task arguments as strings.
    bootstrap = (bootstrap == 'yes')
    # Create a new source distribution as a tarball.
    local('python setup.py sdist --formats=gztar')
    # Figure out the release name and tarball filename.
    dist = local('python setup.py --fullname', capture=True).strip()
    tarball = '%s.tar.gz' % dist
    # Create a place where we can unzip the source tarball.
    run('mkdir /tmp/mutalyzer')
    # Upload the source tarball to the temporary folder on the server.
    put('dist/%s' % tarball, '/tmp/mutalyzer/%s' % tarball)
    # Go to that directory and unzip it.
    with cd('/tmp/mutalyzer'):
        run('tar xzf %s' % tarball)
        # Go to the tarball's contents and do the installation.
        with cd('/tmp/mutalyzer/%s' % dist):
            if bootstrap:
                # Install dependencies.
                sudo('bash extras/pre-install.sh')
            # Install Mutalyzer.
            sudo('python setup.py install')
            if bootstrap:
                # Configure Mutalyzer.
                sudo('bash extras/post-install.sh')
            else:
                # Restart services.
                sudo('bash extras/post-upgrade.sh')
    # Now that all is set up, delete the folder again.
    # I don't like to 'sudo rm -Rf' but since there where files created by
    # root in this directory, we have to.
    sudo('rm -Rf /tmp/mutalyzer')
