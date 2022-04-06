__version__ = "0.0.1"

# Module utilizes internal functions only available from V3.6+
import sys

major, minor = sys.version_info.major, sys.version_info.minor
if major < 3 or (major > 2 and minor < 6):
    raise SystemError(
        "This module requires python3.6 or newer. You are using {}.{}".format(
            major, minor
        )
    )
else:
    print("Python version check: ok")

# Use old formatting syntax to prevent a SyntaxError from running above version check
print("importing pyOER v{} from {}".format(__version__, __file__))

from .elog import *
from .measurement import *
from .calibration import *
from .icpms import *
from .experiment import *
from .sample import *
from .tof import *
from .results_collections import *
from .modelling import *
from .iss import ISS
