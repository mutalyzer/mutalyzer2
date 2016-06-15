"""
Module for retrieving reference files from either the cache or external
services.

Reference gene models and metadata are completely cached in the database. The
sequences themselves are cached on the filesystem.
"""


from __future__ import unicode_literals

from mutalyzer.db.models import Reference
from mutalyzer.parsers import lrg


def load_reference(accession):
    try:
        reference = Reference.query.filter_by(accession=accession).one()
    except NotFound:
        reference = lrg.load_reference(fetch(accession))

    return create_record(reference)


'''
    # Make a filename based upon the identifier.
    filename = self._name_to_file(identifier)

    if not os.path.isfile(filename):
        # We can't find the file.
        filename = self.fetch(identifier)

    if filename is None:
        # Notify batch to skip all instance of identifier.
        self._output.addOutput('BatchFlags', ('S1', identifier))
        # return None in case of error.
        return None

    # Now we have the file, so we can parse it.
    file_handle = bz2.BZ2File(filename, 'r')

    # Create GenRecord.Record from LRG file.
    record = lrg.create_record(file_handle.read())
    file_handle.close()

    # We don't create LRGs from other sources, so id is always the same
    # as source_id.
    record.id = identifier
    record.source_id = identifier

    return record

    def fetch(self, name):
        """
        Fetch the LRG file and store in the cache directory. First try to
        grab the file from the confirmed section, if this fails, get it
        from the pending section.

        :arg unicode name: The name of the LRG file to fetch.

        :returns: the full path to the file; None in case of an error.
        :rtype: unicode
        """
        prefix = settings.LRG_PREFIX_URL
        url = prefix + '{}.xml'.format(name)
        pending_url = prefix + 'pending/{}.xml'.format(name)

        try:
            return self.downloadrecord(url, name)
        except urllib2.URLError:
            # Catch error: file not found.
            pass

        try:
            # Try to get the file from the pending section.
            filename = self.downloadrecord(pending_url, name)
            self._output.addMessage(
                __file__, 2, 'WPEND',
                'Warning: LRG file {} is a pending entry.'.format(name))
            return filename
        except urllib2.URLError:
            self._output.addMessage(
                __file__, 4, 'ERETR', 'Could not retrieve {}.'.format(name))
            # Explicit return in case of an Error.
            return None

    def downloadrecord(self, url, name=None):
        """
        Download an LRG record from an URL.

        :arg unicode url: Location of the LRG record.

        :returns: The full path to the file or Nonein case of failure.
        :rtype: unicode
        """
        lrg_id = name or os.path.splitext(os.path.split(url)[1])[0]
        # if not lrg_id.startswith('LRG'):
        #     return None
        filename = self._name_to_file(lrg_id)

        # TODO: Properly read the file contents to a unicode string and write
        # it utf-8 encoded.
        handle = urllib2.urlopen(url)
        info = handle.info()

        if (info['Content-Type'] == 'application/xml' and
                'Content-length' in info):
            # Looks like a valid LRG file.

            length = int(info['Content-Length'])
            if 512 < length < settings.MAX_FILE_SIZE:
                raw_data = handle.read()
                handle.close()

                # Do an md5 check.
                md5sum = self._calculate_hash(raw_data)
                try:
                    reference = Reference.query.filter_by(
                        accession=lrg_id).one()
                    md5_db = reference.checksum
                except NoResultFound:
                    md5_db = None

                if md5_db is None:
                    # Note: The abstraction seems a bit off here, but we
                    # prefer to set `Reference.source` to `lrg` and not to
                    # `url`, since the former is more specific.
                    reference = Reference(lrg_id, md5sum, 'lrg')
                    session.add(reference)
                    session.commit()
                elif md5sum != md5_db:
                    # Hash has changed for the LRG ID.
                    self._output.addMessage(
                        __file__, -1, 'WHASH',
                        'Warning: Hash of {} changed from {} to {}.'.format(
                            lrg_id, md5_db, md5sum))
                    Reference.query.filter_by(accession=lrg_id).update(
                        {'checksum': md5sum})
                    session.commit()
                else:
                    # Hash the same as in db.
                    pass

                if not os.path.isfile(filename):
                    return self.write(raw_data, lrg_id)
                else:
                    # This can only occur if synchronus calls to mutalyzer are
                    # made to recover a file that did not exist. Still leaves
                    # a window in between the check and the write.
                    return filename
            else:
                self._output.addMessage(
                    __file__, 4, 'EFILESIZE',
                    'Filesize is not within the allowed boundaries.')
        else:
            self._output.addMessage(
                __file__, 4, 'ERECPARSE', 'This is not an LRG record.')
        handle.close()

    def write(self, raw_data, filename):
        """
        Write raw LRG data to a file. The data is parsed before writing,
        if a parse error occurs None is returned.

        :arg str raw_data: The data.
        :arg unicode filename: The intended name of the file.

        :returns: The full path and name of the file written, None in case of
          an error.
        :rtype: unicode
        """
        # Dirty way to test if a file is valid,
        # Parse the file to see if it's a real LRG file.
        try:
            lrg.create_record(raw_data)
        except DOMException:
            self._output.addMessage(
                __file__, 4, 'ERECPARSE', 'Could not parse file.')
            # Explicit return on Error.
            return None

        # Returns full path.
        return self._write(raw_data, filename)
'''
