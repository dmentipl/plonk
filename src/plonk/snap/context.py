"""Context manager for Snaps."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .snap import Snap


class SnapCacheContext:
    """Context manager for caching particle arrays on a Snap.

    Parameters
    ----------
    snap
        The Snap object.
    cache
        Turn caching arrays on or off during a context.
    """

    def __init__(self, snap: Snap, cache: bool):
        self.snap = snap
        self.cache = cache
        self.previous_cache_arrays = snap.cache_arrays
        self.loaded_arrays = snap.loaded_arrays()

    def __enter__(self):
        """Enter context."""
        self.snap.cache_arrays = self.cache

    def __exit__(self, *args):
        """Exit context."""
        self.snap.cache_arrays = self.previous_cache_arrays
        if self.previous_cache_arrays is False:
            self.snap.bulk_unload(
                list(set(self.snap.loaded_arrays()).difference(self.loaded_arrays))
            )


def context(snap: Snap, cache: bool = True):
    """Context manager for caching particle arrays on a Snap.

    Caching arrays in memory improves performance by saving slow
    reads from disc. Caching arrays is turned on by default on a Snap.
    However, it can be useful to turn it off for low memory
    requirements, such as reading from many snapshots. In this case it
    may be useful to use this context manager to improve performance
    inside the context without then leaving those arrays in memory
    afterwards.

    Parameters
    ----------
    cache
        Turn caching arrays on or off during a context.

    Examples
    --------
    Temporarily turn on caching arrays.

    # Caching off
    >>> snap.cache_arrays
    False

    # No loaded arrays
    >>> snap.loaded_arrays()
    ()

    # Produce an image which loads arrays
    >>> with snap.context(cache=True):
    ...     ax = snap.image('density')

    # No loaded arrays
    >>> snap.loaded_arrays()
    ()
    """
    return SnapCacheContext(snap=snap, cache=cache)
