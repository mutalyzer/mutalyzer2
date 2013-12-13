"""
WSGI interface to the Mutalyzer website.

The WSGI interface is exposed through the module variable 'application'.
Static files are not handled by this interface and should be served through
the '/static' url prefix separately.

Example Apache/mod_wsgi configuration:

    Alias /static /var/www/mutalyzer/static
    WSGIScriptAlias / /usr/local/bin/mutalyzer-website

Another common practice is to use Nginx to directly serve the static files
and act as a reverse proxy server to the Mutalyzer HTTP server.

Example Nginx configuration:

    server {
      listen 80;
      location /static/ {
        root /var/www/mutalyzer/static;
        if (-f $request_filename) {
          rewrite ^/static/(.*)$  /static/$1 break;
        }
      }
      location / {
        proxy_read_timeout 300;  # 5 minutes
        proxy_pass http://127.0.0.1:8080;
      }
    }

You can also use the built-in HTTP server by running this file directly. This
will give you a single-threaded server suitable for development which will
also serve the static files.
"""


import argparse
import os
import pkg_resources

from .. import website


application = website.app.wsgifunc()


def debugserver():
    """
    Run the website with the Python built-in HTTP server.
    """
    # There's really no sane way to make web.py serve static files other than
    # providing it with a `static` directory, so we just jump to the template
    # directory where it can find this.
    os.chdir(pkg_resources.resource_filename('mutalyzer', 'templates'))

    website.app.run()


def main():
    """
    Command-line interface to the website.
    """
    parser = argparse.ArgumentParser(
        description='Mutalyzer website.')
    parser.add_argument(
        'port', metavar='NUMBER', type=int, nargs='?', default=8080,
        help='port to run the website on (default: 8080)')

    args = parser.parse_args()
    debugserver()


if __name__ == '__main__':
    main()
