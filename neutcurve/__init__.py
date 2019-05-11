"""
============
neutcurve
============

Fit and draw neutralization curves.

Importing this package imports the following main functions
and classes into the package namespace:

 - :mod:`neutcurve.curvefits.CurveFits`: class to fit curves to
   data, then access and plot results.

 - :class:`neutcurve.hillcurve.HillCurve`: class to fit a single Hill curve.
"""

__author__ = 'Jesse Bloom'
__email__ = 'jbloom@fredhutch.org'
__version__ = '0.2.2'
__url__ = 'https://github.com/jbloomlab/neutcurve'

from neutcurve.curvefits import CurveFits  # noqa: F401
from neutcurve.hillcurve import HillCurve  # noqa: F401
