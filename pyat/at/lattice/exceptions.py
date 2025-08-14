"""Custom Warning and Error subclasses unique for AT.
"""


class AtError(Exception):
    pass


class AtWarning(UserWarning):
    pass
