"""Internal implementation package for the three numbered pipeline commands.

Users should run the numbered files in ``scripts/`` rather than this package directly.
Only the central path object is re-exported here as part of the package's public API.
"""

from .config import PipelinePaths

# ``__all__`` documents what ``from figure1d_pipeline import *`` would expose.
__all__ = ["PipelinePaths"]
