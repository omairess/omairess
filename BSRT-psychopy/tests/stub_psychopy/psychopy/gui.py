"""A dialog that answers with whatever the test queued."""

QUEUE = []


class DlgFromDict(object):
    def __init__(self, dictionary, title='', order=(), tip=None, **kw):
        self.dictionary = dictionary
        if QUEUE:
            dictionary.update(QUEUE.pop(0))
            self.OK = True
        else:
            self.OK = False
