"""Longpipe command-line interface."""

from importlib.metadata import PackageNotFoundError, version


try:
    __version__ = version("longpipe")
except PackageNotFoundError:
    __version__ = "0+unknown"
