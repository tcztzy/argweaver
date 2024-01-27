from importlib.metadata import version

from . import s
from .s import *  # noqa: F403

__doc__ = s.__doc__
__version__ = version("argweavers")

if hasattr(s, "__all__"):
    __all__ = s.__all__
