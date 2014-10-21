String representations
======================

We live in a global economy with many different languages and alphabets. Using
byte strings for text and just assuming everything is ASCII encoded is
suboptimal and *will* lead to bugs. These bugs may even be security issues.

That's why Mutalyzer uses unicode strings wherever possible and tries to be
aware of encodings when dealing with input and output. Here we describe how we
do it.


String representations in Python
--------------------------------

Since Mutalyzer only runs on Python 2.7, we can ignore all older Python versions
and Python 3. So, the two main string types in Python are:

1. `str`, byte strings
2. `unicode`, unicode strings

Byte strings are the default string type in Python 2.7 and are for example the
type you get when writing a string literal::

    >>> type('mutalyzer')
    <type 'str'>

Unicode string literals can be written using the ``u`` prefix::

    >>> type(u'mutalyzer')
    <type 'unicode'>

Many modules from the Python standard library and also third party libraries
consume and produce byte strings by default and may or may not work correctly
with unicode strings.


Unicode strategy
----------------

Internally, all strings should be represented by unicode strings as much as
possible. The main exceptions are large reference sequence strings. These can
often better be BioPython sequence objects, since that is how we usually get
them in the first place.

Our main strategy is as follows:

1. We use ``from __future__ import unicode_literals`` at the top of every
   file.
2. All incoming strings are decoded to unicode (if necessary) as soon as
   possible.
3. Outgoing strings are encoded to UTF8 (if necessary) as late as possible.
4. BioPython sequence objects can be based on byte strings as well as unicode
   strings.
5. In the database, everything is UTF8.
6. We must be aware of the encoding of files supplied by the user or
   downloaded from external sources.

Point 1 ensures that `all string literals in our source code will be unicode
strings <http://python-future.org/unicode_literals.html>`_::

    >>> from __future__ import unicode_literals
    >>> type('mutalyzer')
    <type 'unicode'>

As for point 4, sometimes this may even change under our eyes (e.g., calling
``.reverse_complement()`` will change it to a byte string). We don't care as
long as they're BioPython objects, only when we get the sequence out we must
have it as unicode string. Their contents are always in the ASCII range
anyway.

Although `Bio.Seq.reverse_complement` works fine on Python byte strings (and
we used to rely on that), it crashes on a Python unicode string. So we take
care to only use it on BioPython sequence objects and wrote our own reverse
complement function for unicode strings
(`mutalyzer.util.reverse_complement`).


Files
-----

The Python builtin `open
<https://docs.python.org/2/library/functions.html#open>`_ cannot decode file
contents and just yields byte strings. Therefore, we typically use `io.open
<https://docs.python.org/2/library/io.html#io.open>`_ instead, which accepts
an `encoding` argument.

Downloaded reference files are stored UTF8 encoded (and then bzipped). We can
assume UTF8 encoding when reading them back from disk.

We try to detect the encoding of user uploaded text files (batch jobs, GenBank
files) and assume UTF8 if detection fails.


Libraries
---------

SQLAlchemy, our database toolkit, transparently sends both byte strings and
unicode strings UTF8 encoded to the database and presents all strings as
unicode strings to us.

The webframework Mutalyzer uses, Flask, is also fully `unicode based
<http://flask.pocoo.org/docs/0.10/unicode/>`_.

The Mutalyzer webservices are based on Spyne. The Spyne documentation `has the
following to say <http://spyne.io/docs/2.10/manual/03_types.html#strings>`_
about its `String` and `Unicode` types:

    There are two string types in Spyne: `spyne.model.primitive.Unicode` and
    `spyne.model.primitive.String` whose native types are `unicode` and `str`
    respectively.

    Unlike the Python `str`, the Spyne `String` is not for arbitrary byte
    streams. You should not use it unless you are absolutely, positively sure
    that you need to deal with text data with an unknown encoding. In all
    other cases, you should just use the `Unicode` type. They actually look
    the same from outside, this distinction is made just to properly deal with
    the quirks surrounding Python-2's `unicode` type.

    Remember that you have the `ByteArray` and `File` types at your disposal
    when you need to deal with arbitrary byte streams.

    The `String` type will be just an alias for `Unicode` once Spyne gets
    ported to Python 3. It might even be deprecated and removed in the future,
    so make sure you are using either `Unicode` or `ByteArray` in your
    interface definitions.

So let's not ignore that and never use `String` in our webservice interface.

The pyparsing library is used for parsing HGVS variant descriptions. Overall
it can deal with unicode input and also yields unicode output in that
case. However, there are some exceptions where we explicitely have to decode
to a unicode string (for example, omitted optional parts yield the empty byte
string).


Python 3
--------

The situation in Python 3 is very different from Python 2.7. The two main
string types in Python 3 are:

1. `str`, unicode strings
2. `byte`, byte strings

Unicode strings are the default string type in Python 3 and are for example
the type you get when writing a string literal::

    >>> type('mutalyzer')
    <class 'str'>

Byte string literals can be written using the ``b`` prefix::

    >>> type(b'mutalyzer')
    <class 'bytes'>

Many modules from the Python standard library and also third party libraries
consume and produce unicode strings by default and may or may not work
correctly with byte strings.

What does this mean for Mutalyzer? Actually, our current approach takes us
quite a bit closer to how things are generally done in Python 3. However,
Mutalyzer is very much not Python 3 compatible, even the unicode handling
parts are only valid in Python 2.7 on some points.
