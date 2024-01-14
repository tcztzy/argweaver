from importlib.metadata import version

from . import s

__doc__ = s.__doc__
__version__ = version(__name__)

if hasattr(s, "__all__"):
    __all__ = s.__all__
