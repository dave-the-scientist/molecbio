"""
Several modules to handle running, chaining, and parsing Rosetta tasks.

"""

__version__ = '0.5'
__author__ = 'Dave Curran'
__all__ = ['pdbFilter', '_options']

# # # # #  Version Checks  # # # # #
import sys
version = sys.version_info[:2]
if version < (2, 5):
    raise RuntimeError('This package requires python version 2.5 or newer.')
elif version >= (3, 0):
    import warnings
    warnings.warn('This package was designed for python 2.5 or higher; it may or may not work with version 3.x.', RuntimeWarning)

# # # # #  Package Exceptions and Warnings  # # # # #
class RosettaError(Exception):
    pass
class PathError(RosettaError):
    pass

def OptionWarning(option):
    s = '\tRosetta Option Warning: The %s option could not be found or was not valid. This option is important, and the program may not run unless it is corrected.\n' % option
    sys.stderr.write(s); sys.stderr.flush()
def DirectoryPathWarning(directory, path):
    s = '\tRosetta Directory Path Warning: The %s directory could not be found at %s. This program will not run unless this path is corrected.\n' % (directory, path)
    sys.stderr.write(s); sys.stderr.flush()
def ExecutablePathWarning(executable, path):
    s = '\tRosetta Executable Path Warning: The %s executable could not be found at %s. This program will not run unless this path is corrected.\n' % (executable, path)
    sys.stderr.write(s); sys.stderr.flush()

# # # # #  Package Imports  # # # # #
import _options
options = _options.OptionsDict()
options.save()

import pdbFilter
try: from molecbio.rosettaApps import ScoreView
except ImportError:
    print "Could not import the ScoreView module. This is most likely becuse the 'ttk' dependency could not be found."
from classes import *

