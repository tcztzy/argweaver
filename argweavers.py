from importlib.metadata import version

import argweaver.s

__doc__ = argweaver.s.__doc__
__version__ = version(__file__)

if hasattr(argweaver.s, "__all__"):
    __all__ = argweaver.s.__all__
