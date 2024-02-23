from importlib.metadata import version

from . import argweavers  # make ruff happy
from .argweavers import *  # noqa: F403

__doc__ = argweavers.__doc__
__version__ = version(__package__)

if hasattr(argweavers, "__all__"):
    __all__ = argweavers.__all__
