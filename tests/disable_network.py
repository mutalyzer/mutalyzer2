"""
Disable all non-localhost network access by patching the socket module.

Adapted from Astropy:
https://github.com/astropy/astropy/blob/master/astropy/tests/disable_internet.py
"""


from __future__ import unicode_literals

import contextlib
import socket
import urllib2


# Save original socket methods for restoration. These are global so that
# re-calling the turn_off_network function doesn't overwrite them again.
socket_original = socket.socket
socket_create_connection = socket.create_connection
socket_bind = socket.socket.bind
socket_connect = socket.socket.connect


NETWORK_OFF = False


# urllib2 uses a global variable to cache its default "opener" for opening
# connections for various protocols; we store it off here so we can restore to
# the default after re-enabling network.
_orig_opener = None


def no_network_socket(original_socket):
    """
    Wraps ``original_socket``, which in most cases is assumed to be a
    `socket.socket` method, to raise an `AssertionError` for any operations
    on non-local AF_INET sockets.
    """
    def new_socket(*args, **kwargs):
        if isinstance(args[0], socket.socket):
            if not args[0].family in (socket.AF_INET, socket.AF_INET6):
                # Should be fine in all but some very obscure cases. More to
                # the point, we don't want to affect AF_UNIX sockets.
                return original_socket(*args, **kwargs)
            host = args[1][0]
            addr_arg = 1
            valid_hosts = ('localhost', '127.0.0.1', '::1')
        else:
            # The only other function this is used to wrap currently is
            # socket.create_connection, which should be passed a 2-tuple, but
            # we'll check just in case.
            if not (isinstance(args[0], tuple) and len(args[0]) == 2):
                return original_socket(*args, **kwargs)
            host = args[0][0]
            addr_arg = 0
            valid_hosts = ('localhost', '127.0.0.1')

        hostname = socket.gethostname()
        fqdn = socket.getfqdn()

        if host in (hostname, fqdn):
            host = 'localhost'
            new_addr = (host, args[addr_arg][1])
            args = args[:addr_arg] + (new_addr,) + args[addr_arg + 1:]

        if any([h in host for h in valid_hosts]):
            return original_socket(*args, **kwargs)
        else:
            raise AssertionError('an attempt was made to access the network')

    return new_socket


def turn_off_network(verbose=False):
    """
    Disable network access via Python by preventing connections from being
    created using the socket module.  Presumably this could be worked around by
    using some other means of accessing the internet, but all default python
    modules (urllib, requests, etc.) use socket [citation needed].
    """
    global NETWORK_OFF
    global _orig_opener

    if NETWORK_OFF:
        return

    NETWORK_OFF = True

    __tracebackhide__ = True
    if verbose:
        print('Network access disabled')

    # Update urllib2 to force it not to use any proxies. Must use {} here (the
    # default of None will kick off an automatic search for proxies).
    _orig_opener = urllib2.build_opener()
    no_proxy_handler = urllib2.ProxyHandler({})
    opener = urllib2.build_opener(no_proxy_handler)
    urllib2.install_opener(opener)

    socket.create_connection = no_network_socket(socket_create_connection)
    socket.socket.bind = no_network_socket(socket_bind)
    socket.socket.connect = no_network_socket(socket_connect)

    return socket


def turn_on_network(verbose=False):
    """
    Restore network access. Not used, but kept in case it is needed.
    """
    global NETWORK_OFF
    global _orig_opener

    if not NETWORK_OFF:
        return

    NETWORK_OFF = False

    if verbose:
        print('Network access enabled')

    urllib2.install_opener(_orig_opener)

    socket.create_connection = socket_create_connection
    socket.socket.bind = socket_bind
    socket.socket.connect = socket_connect
    return socket


@contextlib.contextmanager
def no_network(verbose=False):
    """
    Context manager to temporarily disable network access (if not already
    disabled). If it was already disabled before entering the context manager
    (i.e. `turn_off_network` was called previously) then this is a no-op and
    leaves network access disabled until a manual call to `turn_on_network`.
    """
    already_disabled = NETWORK_OFF

    turn_off_network(verbose=verbose)
    try:
        yield
    finally:
        if not already_disabled:
            turn_on_network(verbose=verbose)
