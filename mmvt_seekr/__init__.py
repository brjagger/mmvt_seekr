"""
mmvt_seekr
SEEKR toolkit for multiscale MD/BD/Milestoning simulations
"""

# Add imports here
from .analyze import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
