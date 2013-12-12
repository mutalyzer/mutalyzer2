"""
WSGI interface to the Mutalyzer website.

The WSGI interface is exposed through the module variable 'application'.
Static files are not handled by this interface and should be served through
the '/static' url prefix separately.

Example Apache/mod_wsgi configuration:

    Alias /static /var/www/mutalyzer/static
    WSGIScriptAlias / /usr/local/bin/mutalyzer-website

You can also use the built-in HTTP server by running this file directly.
Note, however, that static files are only found by this server in a 'static'
subdirectory of the current working directory. If you're running Mutalyzer
from its source code directory, you can satisfy this by creating a quick
symbolic link:

    ln -s mutalyzer/templates/static

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
"""


import argparse

from .. import website


application = website.app.wsgifunc()


def debugserver():
    """
    Run the website with the Python built-in HTTP server.
    """
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
