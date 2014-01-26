"""
Simplistic interface to setting a user announcement.

We store up to one announcement at a time, containing a text body and an
optional 'read more' URL. Since getting the announcement should be really
fast, it can be done on every website pageview without problems.
"""


from mutalyzer.redisclient import client


def set_announcement(body, url=None):
    """
    Set announcement.
    """
    announcement = {'body': body}
    if url:
        announcement['url'] = url

    client.hmset('announcement', announcement)


def unset_announcement():
    """
    Unset announcement.
    """
    client.delete('announcement')


def get_announcement():
    """
    Get announcement.
    """
    announcement = client.hgetall('announcement')
    if not announcement:
        return

    return {'body': announcement['body'],
            'url': announcement.get('url')}
