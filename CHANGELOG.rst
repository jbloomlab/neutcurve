=========
Changelog
=========

All notable changes to this project will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com>`_.

2.1.0
-----

Added
+++++
- Added the ``draw_in_bounds`` option to the curve plotting so that curves never are plotted outside the data range. Addresses [this issue](https://github.com/jbloomlab/neutcurve/issues/59).

2.0.1
-----

Fixed
+++++
- Don't throw an uninterpretable error if ``CurveFits.fitParams`` called with no curves, instead just return empty data frame.

2.0.0
-----

Added
+++++
- The curve fitting parameters (top, bottom, slope) can now be constrained to a range in addition to being completely free or fixed. This can help with fitting some curves more sensibly (see [this issue](https://github.com/jbloomlab/neutcurve/issues/53)). Specifically:
  - ``fixtop`` and ``fixbottom`` parameters to ``HillCurve`` can be 2-tuples of bounds
  - added ``fixslope`` parameter to ``HillCurve`` and ``CurveFits``
  - New ``constrain_params_range`` notebook tests and documents this functionality.

- Add ``no_curve_fit_first`` argument to ``HillCurve`` to aid debugging/development.

- Improvements to metrics for assessing curve fit (see [here](https://github.com/jbloomlab/neutcurve/issues/55#issuecomment-2016975219)):
  - The coefficient of determination (``r2``) now is one if all points are fit by a straight line, rather than engative infinity.
  - A root-mean-square-deviation (square root of mean residual) is now calculated as the ``rmsd`` attribute of ``HillCurve`` objects and reported in fit parameter summaries from ``CurveFits``.

1.1.2
-----

Fixed
+++++
- ``CurveFits.plotReplicates`` no longer fails if too many replicates, but instead recycles colors and markers. To make these combinations unique by default, also added another marker to ``CBMARKERS``.

1.1.1
-----

Fixed
+++++
- Fix internal bug in ``HillCurve`` that led to failure to try alternative fitting methods if initial fitting failed as it sometimes does for problematic curves.

1.1.0
-----

Added
+++++
- Report whether curve midpoint is in bounds of data (interpolated) or extrapolated, similar to as was already done for icXX:
  - Added ``HillCurve.midpoint_bound`` and ``HillCurve.midpoint_bound_type`` attributes
  - Report 'midpoint_bound' and 'midpoint_bound_type' in ``CurveFits.fitParams``

1.0.0
-----

Added
+++++
- Changes to improve and change way fitting is done. Main conceptual difference is to first fit curves with fixed slope, then re-fit all parameters:
  - Add ``test_curves.ipynb`` example.
  - Add ``init_slope`` to ``HillCurve`` and ``CurveFits``
  - Add ``fix_slope_first`` to ``HillCurve`` and ``CurveFits`` and make it True by default.

- Add calculation of coefficient of determination for fits (R2) to quantify how well curve fits data:
  - Add ``HillCurve.r2`` attribute.
  - Report coefficient of determination as "r2" in ``CurveFits.fitParams`` output.

- Improve exception handling:
  - Add ``HillCurveFittingError``, and when fitting fails raise this exception rather than ``RuntimeError`` as was done previously.
  - Improve error messages for some exceptions and change a few exception types.

- Add ``allow_reps_unequal_conc`` to ``CurveFits`` to all serum/virus replicates with different concentrations

0.10.1
------

Fixed
+++++
- Fixed bug in ``CurveFits.combineCurveFits`` when using only one fit to combine and specifying one of ``sera``, ``viruses``, or ``serum_virus_replicates_to_drop``

0.10.0
------

Added
+++++
- Added the following parameters to ``CurveFits.combineCurveFits``: ``sera``, ``viruses``, ``serum_virus_replicates_to_drop``.

0.9.0
-----

Added
+++++
- Added the ``CurveFits.combineCurveFits`` static method to combine multiple ``CurveFits`` objects.

0.8.0
-----

Added
+++++
- Added ``attempt_shared_legend`` parameter to ``CurveFits.plotReplicates`` to enable plotting many replicates that differ among sera/viruses.

0.7.0
-----

Added
+++++
- Added ``no_average`` option to ``CurveFits.fitParams``.

Fixed
+++++
- Fixed problem in ``fitParams`` when only some concentrations have replicates.

0.6.0
------

Fixed
+++++
- Code format with ``black``
- Lint with ``ruff``
- Test with GitHub Actions rather than Travis
- Move examples in docs to notebooks so they can be tested with ``nbval`` and added to docs with ``nbsphinx``
- Update minimum Python to 3.8 and test on 3.11

Removed
+++++++
- Stop document the "Rachel-style neutralization" and parsing from Excel as these may be deprecated eventually.

0.5.7
------

Added
+++++
- Return standard errors on fit parameters.

0.5.6
------

Fixed
+++++
- Fixed bug with ``orderlegend`` in ``plotGrid``.

0.5.5
------

Fixed
+++++
- Only import ``dmslogo`` as needed.

0.5.4
-----

Fixed
+++++
- Fixed reading of Excel ``*.xlsx`` files.

0.5.3
-----

Fixed
+++++
- Better fitting of difficult curves by trying multiple optimization methods.

0.5.2
------

Fixed
+++++
- Better error message if virus or serum is `NaN`.

0.5.1
-----

Fixed
++++++
- `CurveFits` now works if `viruses` or `sera` are categorical.

0.5.0
------

Added
++++++
- Added `CurveFits.plotViruses` method.

0.4.2
-----

Fixed
++++++
- Bug fix in ylabel plotting.

0.4.1
------

Fixed
+++++
- Better selection of initial fit parameters when `infectivity_or_neutralized` is 'neutralized'.

0.4.0
------

Added
+++++
- `infectivity_or_neutralized` option to allow fitting of fraction neutralized as well as fraction infectivity.

Fixed
+++++
- `scipy` deprecation warnings.

0.3.1
------

Fixed
++++++
- Fixed bug when IC50 is at lower bound.

0.3.0
-----

Added
+++++
- Ability to draw vertical lines on neutralization curves (`vlines` option to `CurveFits.plotGrid` and `CurveFits.plotSera`).

0.2.5
-----

Fixed
+++++
- Better fit curves that never reach IC50.

0.2.4
-------

Fixed
+++++++
- Fix bug in ymax on some plots generated by `CurveFits`.

0.2.3
-------

Fixed
++++++++
- Fix bug in `CurveFits.plotGrid` when plotting just wildtype.

0.2.2
---------

Added
+++++++
- `ignore_serum_virus` to `CurveFits.plotSera`.

- Added options to `CurveFits.plotGrid` to **not** share x- and y-axis, and to allow different labels.

0.2.1
-------

Added
++++++
- Custom titles for `CurveFits.plotSera`.

0.2.0
-----------

Added
++++++
- Allow exclusion of specific dilutions from *RachelStyle2019* neutralization assays.

- More / better coloring options for `CurveFits.plotSera`.

- Allow more precise sizing of `CurveFits` plots.

Changed
++++++++
- Smaller tick mark sizes.

0.1.0
---------------------------
Initial release

