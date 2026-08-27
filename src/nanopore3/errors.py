"""One base class for every error this package raises at a user.

The command line caught a hand-written tuple of exception types, so any error not
on that list reached the user as a traceback. Two did: picking a consensus backend
whose native tools are not installed, and passing an unknown stage name to `rerun`.
Both had excellent messages and both were presented as a crash.

Every error meant for a person now inherits from :class:`Nanopore3Error` as well as
its natural base, so ``except ValueError`` keeps working where it already did and
the command line has one thing to catch. Anything that does *not* inherit from it
is a bug in this package rather than a mistake by its user, and should reach the
user as a traceback.
"""

from __future__ import annotations


class Nanopore3Error(Exception):
    """Base for every error caused by input, configuration or environment."""


__all__ = ["Nanopore3Error"]
