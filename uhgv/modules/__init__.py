from uhgv.modules import (
    download,
    classify,
)

try:
    from importlib import metadata
except ImportError:
    import importlib_metadata as metadata

__version__ = metadata.version("uhgv")
