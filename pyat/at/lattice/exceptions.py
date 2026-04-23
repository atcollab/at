"""Custom Warning and Error subclasses unique for AT."""

import warnings


class AtError(Exception):
    """AT-specific errors."""
    pass


class AtWarning(UserWarning):
    """AT-specific warnings."""
    pass


warnings.filterwarnings("always", category=AtWarning)
