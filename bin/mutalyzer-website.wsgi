#!/usr/bin/env python

"""
WSGI interface to the Mutalyzer website.

The WSGI interface is exposed through the module variable 'application'.
Static files are not handled by this interface and should be served through
the '/base' url prefix separately.

Example Apache/mod_wsgi configuration:

  Alias /base /var/www/mutalyzer/base
  WSGIScriptAlias / /usr/local/bin/mutalyzer-website.wsgi

You can also use the built-in HTTP server by running this file directly.
Note, however, that static files are not served by this server. A common
pattern is to use Nginx as a proxy and static file server.

Start the built-in HTTP server on port 8080:

  /usr/local/bin/mutalyzer-website.wsgi 8080

Example Nginx configuration:

  server {
    listen 80;
    location /base/ {
      root /var/www/mutalyzer/base;
      if (-f $request_filename) {
        rewrite ^/base/(.*)$  /base/$1 break;
      }
    }
    location / {
      proxy_read_timeout 300;  # 5 minutes
      proxy_pass http://127.0.0.1:8080;
    }
  }

@todo: Integrate webservice.py (http://webpy.org/cookbook/webservice/).
@todo: Move /templates/base to /static for web.py compatibility.
"""


from mutalyzer import website


application = website.app.wsgifunc()


if __name__ == '__main__':
    website.app.run()
