"""
==================
hillcurve
==================

Defines :class:`HillCurve` for fitting neutralization curves.
"""

import collections
import math
import warnings

import matplotlib.pyplot as plt

import numpy

import pandas as pd

import scipy.optimize


class HillCurveFittingError(Exception):
    """Error fitting a :class:`HillCurve`."""

    pass


class HillCurve:
    r"""A fitted Hill curve, optionally with free baselines.

    Fits :math:`f\left(c\right) = b + \frac{t - b}{1 + \left(c/m\right)^s}`
    where :math:`f\left(c\right)` is the fraction infectivity remaining
    at concentration :math:`c`, :math:`m` is the midpoint
    of the neutralization curve, :math:`t` is the top
    value (e.g., 1), :math:`b` is the bottom value (e.g., 0),
    and :math:`s` is the slope of the curve. Because
    :math:`f\left(c\right)` is the fraction infectivity remaining, we expect
    :math:`f\left(c\right)` to get smaller as :math:`c` gets larger.
    This should lead us to fit :math:`s > 0`.

    When :math:`t = 1` and :math:`b = 0`, this equation is identical to the
    `Hill curve <https://en.wikipedia.org/wiki/Hill_equation_(biochemistry)>`_,
    except that we are calculating the fraction **unbound** rather than
    the fraction bound.

    You may want to fit the fraction neutralized rather than the fraction
    infectivity. In that case, set `infectivity_or_neutralized='neutralized'
    and then the equation that is fit will
    be :math:`f\left(c\right) = t + \frac{b - t}{1 + \left(c/m\right)^s}`,
    which means that `f\left(c\right)` gets larger rather than smaller
    as :math:`c` increases.

    Args:
        `cs` (array-like)
            Concentrations of antibody / serum.
        `fs` (array-like)
            Fraction infectivity remaining at each concentration.
        `infectivity_or_neutralized` ({'infectivity', 'neutralized'})
            Fit the fraction infectivity (:math:`f\left(c\right)` decreases as
            :math:`c` increases) or neutralized (:math:`f\left(c\right)`
            increases as :math:`c` increases). See equations above.
        `fix_slope_first` (bool)
            If True, initially fit with fixed slope, and then start from those values
            to re-fit all parameters including slope.
        `fs_stderr` (`None` or array-like)
            If not `None`, standard errors on `fs`.
        `fixbottom` (`False`, float, or tuple/list of length 2)
            If `False`, fit bottom of curves as free parameter. If float, fix bottom of
            curve to that value. If length-2 array, constrain bottom to be in that range.
        `fixtop` (`False`, float, or tuple/list of length 2)
            If `False`, fit top of curves as free parameter. If float, fix top of
            curve to that value. If length-2 array, constrain top to be in that range.
        `fixslope` (`False`, float, or tuple/list of length 2)
            If `False`, fit slope of curves as free parameter. If float, fix slope to
            that value. If length-2 array, constrain slope to be in that range.
        `fitlogc` (bool)
            Do we do the actual fitting on the concentrations or log
            concentrations? Gives equivalent results in principle, but
            fitting to log concentrations may be more efficient in pratice.
        `use_stderr_for_fit` (bool)
            Do we use `fs_stderr` for the fitting, or just for plotting?
            Usually it is a good idea to set to `False` and **not** use
            for fitting if you only have a few replicates, and the standard
            error is often not that accurate and so will weight some
            points much more than others in a way that may not be
            justified.
        `init_slope` (float)
            Initial value of slope used in fitting. If `fixslope` is set to a single
            value, it overrides this parameter. If `fixslope` is set to a range,
            `init_slope` is adjusted up/down to be in that range.
        `no_curve_fit_first` (bool)
            Normally the method first tries to do the optimization with `curve_fit` from
            `scipy`, and if that fails then tries `scipy.optimize.minimize`. If you set
            this to `True`, it skips trying with `curve_fit`. In general, this option is
            for debugging by knowledgeable developers, and you should not be using it
            otherwise.

    Attributes:
        `cs` (numpy array)
            Concentrations, sorted from low to high.
        `fs` (numpy array)
            Fraction infectivity, ordered to match sorted concentrations.
        `fs_stderr` (numpy array or `None`)
            Standard errors on `fs`.
        `bottom` (float)
            Bottom of curve, :math:`b` in equation above.
        `top` (float)
            Top of curve, :math:`t` in equation above.
        `midpoint` (float)
            Midpoint of curve, :math:`m` in equation above. Note
            that the midpoint may **not** be the same as the :meth:`ic50`
            if :math:`t \ne 1` or :math:`b \ne 0`.
        `midpoint_bound` (float)
            Midpoint if it falls within fitted concentrations, otherwise
            the lowest concentration if below that or upper of above it.
        `midpoint_bound_type` (str)
            Is midpoint 'interpolated', or an 'upper' or 'lower' bound.
        `slope` (float)
            Hill slope of curve, :math:`s` in equation above.
        `r2` (float)
            Coefficient of determination indicating how well the curve fits the
            data (https://en.wikipedia.org/wiki/Coefficient_of_determination).
        `rmsd` (float)
            Root mean square deviation of fitted to actual values (square root of mean
            residual).
        `params_stdev` (dict or `None`)
            If standard deviations can be estimated on the fit
            parameters, keyed by 'bottom', 'top', 'midpoint',
            and 'slope' and gives standard erorr on each. Note if
            you have replicates we recommend fitting those separately
            and taking standard error rather than using fit stdev.

    Use the :meth:`ic50` method to get the fitted IC50.
    You can use :meth:`ic50_stdev` to get the estimated standard
    deviation on the IC50, although if you have multiple replicates
    you may be better off just fitting to each separately and
    then taking standard error of the individual IC50s.

    Here are some examples. First, we import the necessary modules:

    .. nbplot::

        >>> import numpy
        >>>
        >>> from neutcurve import HillCurve
        >>> from neutcurve.colorschemes import CBPALETTE

    Now simulate some data:

    .. nbplot::

        >>> m = 0.03
        >>> s = 1.9
        >>> b = 0.1
        >>> t = 1.0
        >>> cs = [0.002 * 2**x for x in range(9)]
        >>> fs = [HillCurve.evaluate(c, m, s, b, t) for c in cs]

    Now fit to these data, and confirm that the fitted values
    are close to the ones used for the simulation:

    .. nbplot::

        >>> neut = HillCurve(cs, fs, fixbottom=False)
        >>> numpy.allclose(neut.midpoint, m, atol=1e-4)
        True
        >>> neut.midpoint_bound == neut.midpoint
        True
        >>> neut.midpoint_bound_type
        'interpolated'
        >>> numpy.allclose(neut.slope, s, atol=1e-3)
        True
        >>> numpy.allclose(neut.top, t, atol=1e-4)
        True
        >>> numpy.allclose(neut.bottom, b, atol=1e-3)
        True
        >>> for key, val in neut.params_stdev.items():
        ...     print(f"{key} = {val:.2g}")
        midpoint = 0.062
        slope = 6.3
        top = 0
        bottom = 0.73
        >>> print(f"IC50: {neut.ic50():.3f} +/- {neut.ic50_stdev():.3f}")
        IC50: 0.034 +/- 0.070

    Since we fit the curve to simulated data where the bottom was
    0.1 rather than 0, the midpoint and IC50 are different. Specifically,
    the IC50 is larger than the midpoint as you have to go past the midpoint
    to get down to 0.5 fraction infectivity:

    .. nbplot::

        >>> neut.ic50() > neut.midpoint
        True
        >>> numpy.allclose(neut.ic50(), 0.0337, atol=1e-4)
        True
        >>> numpy.allclose(0.5, neut.fracinfectivity(neut.ic50()))
        True
        >>> neut.fracinfectivity(neut.midpoint) > 0.5
        True

    Now here is an example where we constrain both the top
    and the bottom (to 1 and 0, respectively) and fit
    the curve. Now the midpoint and IC50 are the same:

    .. nbplot::

        >>> b2 = 0
        >>> t2 = 1
        >>> fs2 = [HillCurve.evaluate(c, m, s, b2, t2) for c in cs]
        >>> neut2 = HillCurve(cs, fs2)
        >>> numpy.allclose(neut2.midpoint, m)
        True
        >>> numpy.allclose(neut2.ic50(), m)
        True

    Now let's fit to concentrations that are all **less**
    than the midpoint, so that we never get to complete neutralization.
    The estimated IC50 is unreliable, and so will be returned as `None` unless
    we call :meth:`ic50` with `method` set to 'bound':

    .. nbplot::

        >>> cs3 = [1e-4 * 2**x for x in range(7)]
        >>> (cs3[-1] < m)
        True
        >>> fs3 = [HillCurve.evaluate(c, m, s, b2, t2) for c in cs3]
        >>> neut3 = HillCurve(cs3, fs3)
        >>> neut3.ic50() is None
        True
        >>> numpy.allclose(neut3.ic50(method='bound'), cs3[-1])
        True
        >>> neut3.midpoint_bound == cs3[-1]
        True
        >>> neut3.midpoint_bound_type
        'lower'

    Note that we can determine if the IC50 is interpolated or an upper
    or lower bound using :meth:`ic50_bound`, and get a nice string
    using :meth:`ic50_str`:

    .. nbplot::

        >>> neut.ic50_bound()
        'interpolated'
        >>> neut3.ic50_bound()
        'lower'
        >>> neut.ic50_str()
        '0.0337'
        >>> neut3.ic50_str()
        '>0.0064'

    We can use the :meth:`dataframe` method to get the measured
    data and fit data at selected points. First, we do this
    just at measured points:

    .. nbplot::

        >>> neut.dataframe('measured').round(3)
           concentration  measurement    fit  stderr
        0          0.002        0.995  0.995     NaN
        1          0.004        0.981  0.981     NaN
        2          0.008        0.932  0.932     NaN
        3          0.016        0.791  0.791     NaN
        4          0.032        0.522  0.522     NaN
        5          0.064        0.272  0.272     NaN
        6          0.128        0.154  0.154     NaN
        7          0.256        0.115  0.115     NaN
        8          0.512        0.104  0.104     NaN

    Then we add in one more point:

    .. nbplot::

        >>> neut.dataframe([0.6]).round(3)
           concentration  measurement    fit  stderr
        0          0.002        0.995  0.995     NaN
        1          0.004        0.981  0.981     NaN
        2          0.008        0.932  0.932     NaN
        3          0.016        0.791  0.791     NaN
        4          0.032        0.522  0.522     NaN
        5          0.064        0.272  0.272     NaN
        6          0.128        0.154  0.154     NaN
        7          0.256        0.115  0.115     NaN
        8          0.512        0.104  0.104     NaN
        9          0.600          NaN  0.103     NaN

    In reality, you'd typically just call :meth:`dataframe` with
    the default argument of 'auto' to get a good range to plot.
    This is done if we call the :meth:`HillCurve.plot` method:

    .. nbplot::

        >>> fig, ax = neut.plot()

    Finally, we confirm that we get the same result regardless
    of whether we fit using the concentrations in linear or
    log space:

    .. nbplot::

        >>> neut_linear = HillCurve(cs, fs, fitlogc=False, fixbottom=False)
        >>> all(numpy.allclose(getattr(neut, attr), getattr(neut_linear, attr))
        ...     for attr in ['top', 'bottom', 'slope', 'midpoint'])
        True

    Demonstrate :meth:`HillCurve.icXX`:

    >>> neut.icXX(0.95) is None
    True
    >>> neut.icXX_str(0.95)
    '>0.512'
    >>> neut.icXX_bound(0.95)
    'lower'
    >>> numpy.allclose(neut.icXX(0.95, method='bound'), neut.cs[-1])
    True
    >>> '{:.4f}'.format(neut.icXX(0.8), 4)
    '0.0896'
    >>> numpy.allclose(0.2, neut.fracinfectivity(neut.icXX(0.8)))
    True

    We can quantify the goodness of fit with :attr:`HillCurve.r2`. For these simulated
    data the fit is perfect (coefficient of determination of 1):

    >>> round(neut.r2, 3)
    1.0

    We can also quantify the goodness of fit with :attr:`HillCurve.rmsd`:
    >>> round(neut.rmsd, 3)
    0.0

    Now fit with bounds on the parameters. First, we make bounds cover the true values:

    >>> neut_bounds_cover = HillCurve(
    ...     cs, fs, fixbottom=(0, 0.2), fixtop=(0.9, 1), fixslope=(1, 2),
    ... )
    >>> numpy.allclose(neut_bounds_cover.midpoint, m, atol=1e-4)
    True
    >>> numpy.allclose(neut_bounds_cover.slope, s, atol=1e-3)
    True
    >>> numpy.allclose(neut_bounds_cover.top, t, atol=1e-4)
    True
    >>> numpy.allclose(neut_bounds_cover.bottom, b, atol=1e-3)
    True
    >>> round(neut_bounds_cover.r2, 3)
    1.0
    >>> round(neut_bounds_cover.rmsd, 3)
    0.0

    Next fit with bounds that do not cover the true parameters:
    >>> neut_bounds_nocover = HillCurve(
    ...     cs, fs, fixbottom=(0, 0.05), fixtop=(0.9, 0.95), fixslope=(1, 1.5),
    ... )
    >>> round(neut_bounds_nocover.midpoint, 2)
    0.04
    >>> round(neut_bounds_nocover.slope, 2)
    1.5
    >>> round(neut_bounds_nocover.top, 2)
    0.95
    >>> round(neut_bounds_nocover.bottom, 2)
    0.05
    >>> round(neut_bounds_nocover.r2, 2)
    0.99
    >>> round(neut_bounds_nocover.rmsd, 3)
    0.045

    Now fit with `infectivity_or_neutralized='neutralized'`, which is useful
    when the signal **increases** rather than decreases with increasing
    concentration (as would be the case if measuring fraction bound rather
    than fraction infectivity).

    .. nbplot::

       >>> neut_opp = HillCurve(cs, [1 - f for f in fs],
       ...                      fixtop=False, fixbottom=False,
       ...                      infectivity_or_neutralized='neutralized')
       >>> numpy.allclose(neut_opp.top, 0.9, atol=1e-3)
       True
       >>> numpy.allclose(neut_opp.bottom, 0, atol=1e-3)
       True
       >>> numpy.allclose(neut_opp.midpoint, m)
       True
       >>> neut_opp.ic50() < neut_opp.midpoint
       True
       >>> fig, ax = neut_opp.plot(ylabel='fraction neutralized')

    For internal testing purposes, try with `no_curve_fit_first=False`.

    >>> neut_ncf = HillCurve(cs, fs, fixbottom=False, no_curve_fit_first=True)
    >>> numpy.allclose(neut_ncf.midpoint, m, atol=1e-4)
    True
    >>> numpy.allclose(neut_ncf.slope, s, atol=5e-3)
    True
    >>> numpy.allclose(neut_ncf.top, t, atol=1e-4)
    True
    >>> numpy.allclose(neut_ncf.bottom, b, atol=1e-3)
    True

    >>> neut_bounds_cover_ncf = HillCurve(
    ...     cs,
    ...     fs,
    ...     fixbottom=(0, 0.2),
    ...     fixtop=(0.9, 1),
    ...     fixslope=(1, 2),
    ...     no_curve_fit_first=True,
    ... )
    >>> numpy.allclose(neut_bounds_cover_ncf.midpoint, m, atol=1e-4)
    True
    >>> numpy.allclose(neut_bounds_cover_ncf.slope, s, atol=1e-3)
    True
    >>> numpy.allclose(neut_bounds_cover_ncf.top, t, atol=1e-4)
    True
    >>> numpy.allclose(neut_bounds_cover_ncf.bottom, b, atol=1e-3)
    True
    >>> round(neut_bounds_cover_ncf.r2, 3)
    1.0

    >>> neut_bounds_nocover_ncf = HillCurve(
    ...     cs,
    ...     fs,
    ...     fixbottom=(0, 0.05),
    ...     fixtop=(0.9, 0.95),
    ...     fixslope=(1, 1.5),
    ...     no_curve_fit_first=True,
    ... )
    >>> round(neut_bounds_nocover_ncf.midpoint, 2)
    0.04
    >>> round(neut_bounds_nocover_ncf.slope, 2)
    1.5
    >>> round(neut_bounds_nocover_ncf.top, 2)
    0.95
    >>> round(neut_bounds_nocover_ncf.bottom, 2)
    0.05
    >>> round(neut_bounds_nocover_ncf.r2, 2)
    0.99

    """

    def __init__(
        self,
        cs,
        fs,
        *,
        infectivity_or_neutralized="infectivity",
        fix_slope_first=True,
        fs_stderr=None,
        fixbottom=0,
        fixtop=1,
        fixslope=False,
        fitlogc=False,
        use_stderr_for_fit=False,
        init_slope=1.5,
        no_curve_fit_first=False,
    ):
        """See main class docstring."""
        # get data into arrays sorted by concentration
        self.cs = numpy.array(cs)
        self.fs = numpy.array(fs)
        if fs_stderr is not None:
            self.fs_stderr = numpy.array(fs_stderr)
            self.fs_stderr = self.fs_stderr[self.cs.argsort()]
        else:
            self.fs_stderr = None
        self.fs = self.fs[self.cs.argsort()]
        self.cs = self.cs[self.cs.argsort()]

        if infectivity_or_neutralized == "infectivity":
            if self.fs[0] < self.fs[-1] and (self.fs[0] < 0.3 and self.fs[-1] > 0.7):
                warnings.warn(
                    "`f` increases with concentration, consider "
                    '`infectivity_or_neutralized="neutralized"'
                )
        elif infectivity_or_neutralized == "neutralized":
            if self.fs[0] > self.fs[-1] and (self.fs[0] > 0.7 and self.fs[-1] < 0.3):
                warnings.warn(
                    "`f` decreases with concentration, consider "
                    '`infectivity_or_neutralized="infectivity"'
                )
        else:
            raise ValueError(f"invalid {infectivity_or_neutralized=}")
        self._infectivity_or_neutralized = infectivity_or_neutralized

        if any(self.cs <= 0):
            raise ValueError(f"concentrations must all be > 0\n{self.cs=}")

        # create initial guess of `(midpoint, slope, bottom, top)`
        # make initial guess for top
        if fixtop is False:
            top = max(1, self.fs.max())
        elif isinstance(fixtop, (int, float)):
            fixtop = float(fixtop)
            top = fixtop
        elif hasattr(fixtop, "__len__") and len(fixtop) == 2:
            fixtop = (fixtop[0], fixtop[1])
            if fixtop[1] <= fixtop[0]:
                raise ValueError(f"invalid {fixtop=}: first element must be < second")
            top = max(fixtop[0], min(max(1, self.fs.max()), fixtop[1]))
            assert fixtop[0] <= top <= fixtop[1]
        else:
            raise ValueError(f"{fixtop=} is not `False`, a number, or length-2 array")
        # make initial guess for bottom
        if fixbottom is False:
            bottom = min(0, self.fs.min())
        elif isinstance(fixbottom, (int, float)):
            fixbottom = float(fixbottom)
            bottom = fixbottom
        elif hasattr(fixbottom, "__len__") and len(fixbottom) == 2:
            fixbottom = (fixbottom[0], fixbottom[1])
            if fixbottom[1] <= fixbottom[0]:
                raise ValueError(
                    f"invalid {fixbottom=}: first element must be < second"
                )
            bottom = max(fixbottom[0], min(min(0, self.fs.min()), fixbottom[1]))
            assert fixbottom[0] <= bottom <= fixbottom[1]
        else:
            raise ValueError(
                f"{fixbottom=} is not `False`, a number, or length-2 array"
            )
        # make initial guess for midpoint
        # if midpoint guess outside range, guess outside range by amount
        # equal to spacing of last two points
        midval = (top - bottom) / 2.0
        if (self.fs > midval).all():
            midpoint = {
                "infectivity": self.cs[-1] ** 2 / self.cs[-2],
                "neutralized": self.cs[0] / (self.cs[-1] / self.cs[-2]),
            }[self._infectivity_or_neutralized]
        elif (self.fs <= midval).all():
            midpoint = {
                "neutralized": self.cs[-1] ** 2 / self.cs[-2],
                "infectivity": self.cs[0] / (self.cs[-1] / self.cs[-2]),
            }[self._infectivity_or_neutralized]
        else:
            # get first index where f crosses midpoint
            i = numpy.argmax((self.fs > midval)[:-1] != (self.fs > midval)[1:])
            assert (self.fs[i] > midval) != (self.fs[i + 1] > midval)
            midpoint = (self.cs[i] + self.cs[i + 1]) / 2.0
        # adjust initial slope if inconsistent with `fixslope`
        if fixslope is False:
            slope = float(init_slope)
        elif isinstance(fixslope, (int, float)):
            fixslope = float(fixslope)
            slope = fixslope
        elif hasattr(fixslope, "__len__") and len(fixslope) == 2:
            fixslope = (fixslope[0], fixslope[1])
            if fixslope[1] <= fixslope[0]:
                raise ValueError(f"invalid {fixslope=}: first element must be < second")
            slope = max(fixslope[0], min(init_slope, fixslope[1]))
            assert fixslope[0] <= slope <= fixslope[1]

        init_tup = (midpoint, slope, bottom, top)

        # first try to fit using curve_fit
        try:
            if no_curve_fit_first:
                raise RuntimeError("skipping curve_fit")
            if fix_slope_first:
                fix_first_init_tup, _ = self._fit_curve(
                    fixtop=fixtop,
                    fixbottom=fixbottom,
                    fitlogc=fitlogc,
                    use_stderr_for_fit=use_stderr_for_fit,
                    init_tup=init_tup,
                    fixslope=slope,
                )
            else:
                fix_first_init_tup = init_tup
            fit_tup, self.params_stdev = self._fit_curve(
                fixtop=fixtop,
                fixbottom=fixbottom,
                fitlogc=fitlogc,
                use_stderr_for_fit=use_stderr_for_fit,
                init_tup=fix_first_init_tup,
                fixslope=fixslope,
            )
        # A RuntimeError is raised by scipy if curve_fit fails:
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
        except RuntimeError:
            # curve_fit failed, try using minimize
            for method in ["TNC", "L-BFGS-B", "SLSQP", "Powell"]:
                if fix_slope_first:
                    fix_first_init_tup = self._minimize_fit(
                        fixtop=fixtop,
                        fixbottom=fixbottom,
                        fitlogc=fitlogc,
                        use_stderr_for_fit=use_stderr_for_fit,
                        method=method,
                        init_tup=init_tup,
                        fixslope=slope,
                    )
                    if fix_first_init_tup is False:
                        continue
                else:
                    fix_first_init_tup = init_tup
                fit_tup = self._minimize_fit(
                    fixtop=fixtop,
                    fixbottom=fixbottom,
                    fitlogc=fitlogc,
                    use_stderr_for_fit=use_stderr_for_fit,
                    method=method,
                    init_tup=fix_first_init_tup,
                    fixslope=fixslope,
                )
                self.params_stdev = None  # can't estimate errors
                if fit_tup is not False:
                    break
            else:
                raise HillCurveFittingError(f"fit failed:\n{self.cs=}\n{self.fs=}")

        for i, param in enumerate(["midpoint", "slope", "bottom", "top"]):
            setattr(self, param, float(fit_tup[i]))
        # debugging check
        for name, fixval, val in [
            ("slope", fixslope, self.slope),
            ("top", fixtop, self.top),
            ("bottom", fixbottom, self.bottom),
        ]:
            if not (
                (fixval is False)
                or (isinstance(fixval, float) and numpy.allclose(fixval, val))
                or (isinstance(fixval, tuple) and (fixval[0] <= val <= fixval[1]))
            ):
                raise ValueError(f"fix{name}={fixval}, but final {name}={val}")

        if self.cs[0] <= self.midpoint <= self.cs[-1]:
            self.midpoint_bound = self.midpoint
            self.midpoint_bound_type = "interpolated"
        elif self.midpoint < self.cs[0]:
            self.midpoint_bound = float(self.cs[0])
            self.midpoint_bound_type = "upper"
        else:
            assert self.midpoint > self.cs[-1]
            self.midpoint_bound = float(self.cs[-1])
            self.midpoint_bound_type = "lower"

        # compute coefficient of determination
        # https://en.wikipedia.org/wiki/Coefficient_of_determination
        assert len(self.cs) == len(self.fs)
        sstot = float(((self.fs - numpy.average(self.fs)) ** 2).sum())
        ssres = float(
            (
                (numpy.array([self.fracinfectivity(c) for c in self.cs]) - self.fs) ** 2
            ).sum()
        )
        if sstot == 0:
            if ssres == 0:
                self.r2 = 1.0
            else:
                self.r2 = 0.0
        else:
            self.r2 = 1.0 - ssres / sstot

        # compute rmsd
        self.rmsd = math.sqrt(ssres / len(self.cs))

    def _fit_curve(
        self,
        *,
        fixtop,
        fixbottom,
        fitlogc,
        use_stderr_for_fit,
        init_tup,
        fixslope,
    ):
        """curve_fit, return `(midpoint, slope, bottom, top), params_stdev`."""

        midpoint, slope, bottom, top = init_tup

        # set up function and initial guesses
        if fitlogc:
            evalfunc = self._evaluate_log
            xdata = numpy.log(self.cs)
            midpoint = numpy.log(midpoint)
        else:
            evalfunc = self.evaluate
            xdata = self.cs

        top_bounds = (-numpy.inf, numpy.inf) if (fixtop is False) else fixtop
        bottom_bounds = (-numpy.inf, numpy.inf) if (fixbottom is False) else fixbottom
        slope_bounds = (-numpy.inf, numpy.inf) if (fixslope is False) else fixslope

        if isinstance(fixslope, float):
            assert fixslope == slope
            if (not isinstance(fixtop, float)) and (not isinstance(fixbottom, float)):
                initguess = [midpoint, bottom, top]
                # https://stackoverflow.com/a/8081580
                bounds = tuple(
                    zip(*[(-numpy.inf, numpy.inf), bottom_bounds, top_bounds])
                )

                def func(c, m, b, t):
                    return evalfunc(c, m, slope, b, t, self._infectivity_or_neutralized)

            elif not isinstance(fixtop, float):
                assert isinstance(fixbottom, float)
                initguess = [midpoint, top]
                bounds = tuple(zip(*[(-numpy.inf, numpy.inf), top_bounds]))

                def func(c, m, t):
                    return evalfunc(
                        c,
                        m,
                        slope,
                        bottom,
                        t,
                        self._infectivity_or_neutralized,
                    )

            elif not isinstance(fixbottom, float):
                assert isinstance(fixtop, float)
                initguess = [midpoint, bottom]
                bounds = tuple(zip(*[(-numpy.inf, numpy.inf), bottom_bounds]))

                def func(c, m, b):
                    return evalfunc(
                        c,
                        m,
                        slope,
                        b,
                        top,
                        self._infectivity_or_neutralized,
                    )

            else:
                assert isinstance(fixtop, float) and isinstance(fixbottom, float)
                initguess = [midpoint]
                bounds = ([-numpy.inf], [numpy.inf])

                def func(c, m):
                    return evalfunc(
                        c,
                        m,
                        slope,
                        bottom,
                        top,
                        self._infectivity_or_neutralized,
                    )

        else:
            assert (fixslope is False) or (
                isinstance(fixslope, tuple) and len(fixslope) == 2
            )
            if (not isinstance(fixtop, float)) and (not isinstance(fixbottom, float)):
                initguess = [midpoint, slope, bottom, top]
                bounds = tuple(
                    zip(
                        *[
                            (-numpy.inf, numpy.inf),
                            slope_bounds,
                            bottom_bounds,
                            top_bounds,
                        ]
                    )
                )

                def func(c, m, s, b, t):
                    return evalfunc(c, m, s, b, t, self._infectivity_or_neutralized)

            elif not isinstance(fixtop, float):
                assert isinstance(fixbottom, float)
                initguess = [midpoint, slope, top]
                bounds = tuple(
                    zip(*[(-numpy.inf, numpy.inf), slope_bounds, top_bounds])
                )

                def func(c, m, s, t):
                    return evalfunc(
                        c, m, s, bottom, t, self._infectivity_or_neutralized
                    )

            elif not isinstance(fixbottom, float):
                assert isinstance(fixtop, float)
                initguess = [midpoint, slope, bottom]
                bounds = tuple(
                    zip(
                        *[
                            (-numpy.inf, numpy.inf),
                            slope_bounds,
                            bottom_bounds,
                        ]
                    )
                )

                def func(c, m, s, b):
                    return evalfunc(c, m, s, b, top, self._infectivity_or_neutralized)

            else:
                assert isinstance(fixtop, float) and isinstance(fixbottom, float)
                initguess = [midpoint, slope]
                bounds = tuple(zip(*[(-numpy.inf, numpy.inf), slope_bounds]))

                def func(c, m, s):
                    return evalfunc(
                        c,
                        m,
                        s,
                        bottom,
                        top,
                        self._infectivity_or_neutralized,
                    )

        assert len(initguess) == len(bounds[0]) == len(bounds[1])
        assert all(lb < ub for (lb, ub) in zip(bounds[0], bounds[1])), bounds
        (popt, pcov) = scipy.optimize.curve_fit(
            f=func,
            xdata=xdata,
            ydata=self.fs,
            p0=initguess,
            sigma=self.fs_stderr if use_stderr_for_fit else None,
            absolute_sigma=True,
            bounds=bounds,
            maxfev=1000,
        )

        perr = numpy.sqrt(numpy.diag(pcov))

        if fitlogc:
            midpoint = numpy.exp(midpoint)

        midpoint = popt[0]
        if isinstance(fixslope, float):
            params_stderr = {"midpoint": perr[0], "slope": 0, "top": 0, "bottom": 0}
            if (not isinstance(fixbottom, float)) and (not isinstance(fixtop, float)):
                bottom = popt[1]
                params_stderr["bottom"] = perr[1]
                top = popt[2]
                params_stderr["top"] = perr[2]
            elif not isinstance(fixbottom, float):
                bottom = popt[1]
                params_stderr["bottom"] = perr[1]
            elif not isinstance(fixtop, float):
                top = popt[1]
                params_stderr["top"] = perr[1]
        else:
            slope = popt[1]
            params_stderr = {
                "midpoint": perr[0],
                "slope": perr[1],
                "top": 0,
                "bottom": 0,
            }
            if (not isinstance(fixbottom, float)) and (not isinstance(fixtop, float)):
                bottom = popt[2]
                params_stderr["bottom"] = perr[2]
                top = popt[3]
                params_stderr["top"] = perr[3]
            elif not isinstance(fixbottom, float):
                bottom = popt[2]
                params_stderr["bottom"] = perr[2]
            elif not isinstance(fixtop, float):
                top = popt[2]
                params_stderr["top"] = perr[2]

        return (midpoint, slope, bottom, top), params_stderr

    def _minimize_fit(
        self,
        *,
        fixtop,
        fixbottom,
        fitlogc,
        use_stderr_for_fit,
        method,
        init_tup,
        fixslope,
    ):
        """Return `(midpoint, slope, bottom, top)` if succeeds or `False` if fails."""

        midpoint, slope, bottom, top = init_tup

        # set up function and initial guesses
        if fitlogc:
            evalfunc = self._evaluate_log
            xdata = numpy.log(self.cs)
            midpoint = numpy.log(midpoint)
            bounds = [(None, None)]
        else:
            evalfunc = self.evaluate
            xdata = self.cs
            bounds = [(0, None)]

        top_bounds = (None, None) if (fixtop is False) else fixtop
        bottom_bounds = (None, None) if (fixbottom is False) else fixbottom
        slope_bounds = (None, None) if (fixslope is False) else fixslope

        if isinstance(fixslope, float):
            if (not isinstance(fixtop, float)) and (not isinstance(fixbottom, float)):
                initguess = [midpoint, bottom, top]
                bounds = bounds + [bottom_bounds, top_bounds]

                def func(c, m, b, t):
                    return evalfunc(c, m, slope, b, t, self._infectivity_or_neutralized)

            elif not isinstance(fixtop, float):
                assert isinstance(fixbottom, float)
                initguess = [midpoint, top]
                bounds.append(top_bounds)

                def func(c, m, t):
                    return evalfunc(
                        c,
                        m,
                        slope,
                        bottom,
                        t,
                        self._infectivity_or_neutralized,
                    )

            elif not isinstance(fixbottom, float):
                assert isinstance(fixtop, float)
                initguess = [midpoint, bottom]
                bounds.append(bottom_bounds)

                def func(c, m, b):
                    return evalfunc(
                        c,
                        m,
                        slope,
                        b,
                        top,
                        self._infectivity_or_neutralized,
                    )

            else:
                assert isinstance(fixtop, float) and isinstance(fixbottom, float)
                initguess = [midpoint]

                def func(c, m):
                    return evalfunc(
                        c,
                        m,
                        slope,
                        bottom,
                        top,
                        self._infectivity_or_neutralized,
                    )

        else:
            assert isinstance(fixslope, tuple) or (fixslope is False)
            bounds.append(slope_bounds)

            if (not isinstance(fixtop, float)) and (not isinstance(fixbottom, float)):
                initguess = [midpoint, slope, bottom, top]
                bounds = bounds + [bottom_bounds, top_bounds]

                def func(c, m, s, b, t):
                    return evalfunc(c, m, s, b, t, self._infectivity_or_neutralized)

            elif not isinstance(fixtop, float):
                assert isinstance(fixbottom, float)
                initguess = [midpoint, slope, top]
                bounds.append(top_bounds)

                def func(c, m, s, t):
                    return evalfunc(
                        c, m, s, bottom, t, self._infectivity_or_neutralized
                    )

            elif not isinstance(fixbottom, float):
                assert isinstance(fixtop, float)
                initguess = [midpoint, slope, bottom]
                bounds.append(bottom_bounds)

                def func(c, m, s, b):
                    return evalfunc(c, m, s, b, top, self._infectivity_or_neutralized)

            else:
                assert isinstance(fixtop, float) and isinstance(fixbottom, float)
                initguess = [midpoint, slope]

                def func(c, m, s):
                    return evalfunc(
                        c,
                        m,
                        s,
                        bottom,
                        top,
                        self._infectivity_or_neutralized,
                    )

        def min_func(p):
            """Evaluate to zero when perfect fit."""
            if (use_stderr_for_fit is None) or (self.fs_stderr is None):
                return sum((func(xdata, *p) - self.fs) ** 2)
            else:
                return sum((func(xdata, *p) - self.fs / self.fs_stderr) ** 2)

        initguess = numpy.array(initguess, dtype="float")
        assert len(initguess) == len(bounds)
        res = scipy.optimize.minimize(min_func, initguess, bounds=bounds, method=method)

        if not res.success:
            return False

        if fitlogc:
            midpoint = numpy.exp(midpoint)

        midpoint = res.x[0]
        if isinstance(fixslope, float):
            if (not isinstance(fixbottom, float)) and (not isinstance(fixtop, float)):
                bottom = res.x[1]
                top = res.x[2]
            elif not isinstance(fixbottom, float):
                bottom = res.x[1]
            elif not isinstance(fixtop, float):
                top = res.x[1]
        else:
            assert (fixslope is False) or isinstance(fixslope, tuple)
            slope = res.x[1]
            if (not isinstance(fixbottom, float)) and (not isinstance(fixtop, float)):
                bottom = res.x[2]
                top = res.x[3]
            elif not isinstance(fixbottom, float):
                bottom = res.x[2]
            elif not isinstance(fixtop, float):
                top = res.x[2]

        return (midpoint, slope, bottom, top)

    def icXX(self, fracneut, *, method="interpolate"):
        """Generalizes :meth:`HillCurve.ic50` to arbitrary frac neutralized.

        For instance, set `fracneut` to 0.95 if you want the IC95, the
        concentration where 95% is neutralized.

        Args:
            `fracneut` (float)
                Compute concentration at which `fracneut` of the virus
                is expected to be neutralized. Note that this is the
                expected fraction **neutralized**, not the fraction
                infectivity.
            `method` (str)
                Can have following values:

                  - 'interpolate': only return a number for ICXX if it
                    is in range of concentrations, otherwise return `None`.

                  - 'bound': if ICXX is out of range of concentrations,
                    return upper or lower measured concentration depending
                    on if ICXX is above or below range of concentrations.
                    Assumes infectivity decreases with concentration.

        Returns:
            Number giving ICXX or `None` (depending on value of `method`).

        """
        fracinf = 1 - fracneut
        if self.top < fracinf and self.bottom < fracinf:
            bound = "lower"
        elif self.top >= fracinf and self.bottom >= fracinf:
            bound = "upper"
        else:
            icXX = self.midpoint * ((self.top - fracinf) / (fracinf - self.bottom)) ** (
                1.0 / self.slope
            )
            if self.cs[0] <= icXX <= self.cs[-1]:
                return icXX
            elif icXX < self.cs[0]:
                bound = "lower"
            else:
                bound = "upper"

        if method == "bound":
            if bound == "upper":
                return self.cs[-1]
            elif bound == "lower":
                return self.cs[0]
            else:
                raise ValueError(f"invalid {bound=}")
        elif method == "interpolate":
            return None
        else:
            raise ValueError(f"invalid {method=}")

    def ic50(self, method="interpolate"):
        r"""IC50 value.

        Concentration where infectivity remaining is 0.5. Equals
        `midpoint` if and only if `top = 1` and `bottom = 0`. Calculated
        from :math:`0.5 = b + \frac{t - b}{1 + (ic50/m)^s}`, which solves to
        :math:`ic50 = m \times \left(\frac{t - 0.5}{0.5 - b}\right)^{1/s}`

        Args:
            `method` (str)
                Can have following values:

                  - 'interpolate': only return a number for IC50 if it
                    is in range of concentrations, otherwise return `None`.

                  - 'bound': if IC50 is out of range of concentrations,
                    return upper or lower measured concentration depending
                    on if IC50 is above or below range of concentrations.
                    Assumes infectivity decreases with concentration.

        Returns:
            Number giving IC50 or `None` (depending on value of `method`).

        """
        return self.icXX(0.5, method=method)

    def ic50_stdev(self):
        r"""Get standard deviation of fit IC50 parameter.

        Calculated just from estimated standard deviation on `midpoint`.
        Note if you have replicates, we recommend fitting separately
        and calculating standard error from those fits rather than
        using this value.

        Returns:
            A number giving the standard deviation, or `None` if cannot
            be estimated or if IC50 is at bound.

        """
        ic50 = self.ic50()
        if ic50 is None:
            return None
        midpoint_stdev = self.params_stdev["midpoint"]
        if midpoint_stdev is None:
            return None
        else:
            return midpoint_stdev * ic50 / self.midpoint

    def icXX_bound(self, fracneut):
        """Like :meth:`HillCurve.ic50_bound` for arbitrary frac neutralized."""
        if self.icXX(fracneut, method="interpolate") is not None:
            return "interpolated"
        else:
            icXX = self.icXX(fracneut, method="bound")
            if icXX == self.cs[0]:
                return "upper"
            elif icXX == self.cs[-1]:
                return "lower"
            else:
                raise ValueError(f"icXX not bound for {fracneut=}")

    def ic50_bound(self):
        """Is IC50 'interpolated', or an 'upper' or 'lower' bound."""
        return self.icXX_bound(0.5)

    def icXX_str(self, fracneut, *, precision=3):
        """Like :meth:`HillCurve.ic50_str` for arbitrary frac neutralized."""
        icXX = f"{{:.{precision}g}}".format(self.icXX(fracneut, method="bound"))
        prefix = {"interpolated": "", "upper": "<", "lower": ">"}[
            self.icXX_bound(fracneut)
        ]
        return f"{prefix}{icXX}"

    def ic50_str(self, precision=3):
        """IC50 as string indicating upper / lower bounds with > or <.

        Args:
            Number of significant digits in returned string.

        """
        return self.icXX_str(0.5, precision=precision)

    def fracinfectivity(self, c):
        """Fraction infectivity at `c` for fitted parameters."""
        return self.evaluate(
            c,
            self.midpoint,
            self.slope,
            self.bottom,
            self.top,
            self._infectivity_or_neutralized,
        )

    @staticmethod
    def evaluate(c, m, s, b, t, infectivity_or_neutralized="infectivity"):
        r""":math:`f\left(c\right) = b + \frac{t-b}{1+\left(c/m\right)^s}`.

        If `infectivity_or_neutralized` is 'neutralized' rather than
        'infectivity', instead return
        :math:`f\left(c\right) = t + \frac{b-t}{1+\left(c/m\right)^s}`.

        """
        if infectivity_or_neutralized == "infectivity":
            return b + (t - b) / (1 + (c / m) ** s)
        elif infectivity_or_neutralized == "neutralized":
            return t + (b - t) / (1 + (c / m) ** s)
        else:
            raise ValueError(f"invalid {infectivity_or_neutralized=}")

    @staticmethod
    def _evaluate_log(logc, logm, s, b, t, infectivity_or_neutralized="infectivity"):
        """Like :class:`HillCurve.evaluate` but on log concentration scale."""
        if infectivity_or_neutralized == "infectivity":
            return b + (t - b) / (1 + numpy.exp(s * (logc - logm)))
        elif infectivity_or_neutralized == "neutralized":
            return t + (b - t) / (1 + numpy.exp(s * (logc - logm)))
        else:
            raise ValueError(f"invalid {infectivity_or_neutralized=}")

    def plot(
        self,
        *,
        concentrations="auto",
        ax=None,
        xlabel="concentration",
        ylabel="fraction infectivity",
        color="black",
        marker="o",
        markersize=6,
        linewidth=1,
        linestyle="-",
        yticklocs=None,
        draw_in_bounds=False,
    ):
        """Plot the neutralization curve.

        Args:
            `concentrations`
                Concentrations to plot, same meaning as for
                :meth:`HillCurve.dataframe`.
            `ax` (`None` or matplotlib axes.Axes object)
                Use to plot on an existing axis. If using an existing
                axis, do **not** re-scale the axis limits to the data.
            `xlabel` (str or `None`)
                Label for x-axis.
            `ylabel` (str or `None`)
                Label for y-axis.
            `color` (str)
                Color of line and point.
            `marker` (str)
                Marker shape: https://matplotlib.org/api/markers_api.html
            `markersize` (float)
                Size of point marker.
            `linewidth` (float)
                Width of line.
            `linestyle` (str)
                Line style.
            `yticklocs` (`None` or list)
                Exact locations to place yticks; `None` means auto-locate.
            `draw_in_bounds` (bool)
                By default, the plotted curve extends a bit outside the actual
                data points on both sides. If this is set to `True`, then the
                plotted curve stops at the bounds of the data points.

        Returns:
            The 2-tuple `(fig, ax)` giving the matplotlib figure and axis.

        """
        data = self.dataframe(concentrations, extend=0 if draw_in_bounds else 0.1)

        if ax is None:
            fig, ax = plt.subplots()
            fig.set_size_inches((4, 3))
            check_ybounds = True
            ylowerbound = -0.05
            yupperbound = 1.05
            ax.autoscale(True, "both")
        else:
            fig = ax.get_figure()
            ax.autoscale(False, "both")
            check_ybounds = False

        ax.plot(
            "concentration",
            "fit",
            data=data,
            linestyle=linestyle,
            linewidth=linewidth,
            color=color,
        )

        ax.errorbar(
            x="concentration",
            y="measurement",
            yerr="stderr",
            data=data,
            fmt=marker,
            color=color,
            markersize=markersize,
            capsize=markersize / 1.5,
        )

        ax.set_xscale("log")
        ax.set_xlabel(xlabel, fontsize=15)
        ax.set_ylabel(ylabel, fontsize=15)
        ax.tick_params("both", labelsize=12, length=5, width=1)
        ax.minorticks_off()
        if yticklocs is not None:
            ax.set_yticks(yticklocs)

        if check_ybounds:
            ymin, ymax = ax.get_ylim()
            ax.set_ylim(min(ymin, ylowerbound), max(ymax, yupperbound))

        return fig, ax

    def dataframe(self, concentrations="auto", extend=0.1):
        """Get data frame with curve data for plotting.

        Useful if you want to get both the points and the fit
        curve to plot.

        Args:
            `concentrations` (array-like or 'auto' or 'measured')
                Concentrations for which we compute the fit values.
                If 'auto' the automatically computed from data
                range using :func:`concentrationRange`. If
                'measured' then only include measured values.
            `extend` (float)
                The value passed to :func:`concentrationRange` if `concentrations`
                is 'auto'.

        Returns:
            A pandas DataFrame with the following columns:

              - 'concentration': concentration
              - 'fit': curve fit value at this point
              - 'measurement': value of measurement at this point,
                or numpy.nan if no measurement here.
              - 'stderr': standard error of measurement if provided.

        """
        if concentrations == "auto":
            concentrations = concentrationRange(self.cs[0], self.cs[-1], extend=extend)
        elif concentrations == "measured":
            concentrations = []
        concentrations = numpy.concatenate([self.cs, concentrations])
        n = len(concentrations)

        points = numpy.concatenate([self.fs, numpy.full(n - len(self.fs), numpy.nan)])

        if self.fs_stderr is None:
            stderr = numpy.full(n, numpy.nan)
        else:
            stderr = numpy.concatenate(
                [self.fs_stderr, numpy.full(n - len(self.fs), numpy.nan)]
            )

        fit = numpy.array([self.fracinfectivity(c) for c in concentrations])

        return (
            pd.DataFrame.from_dict(
                collections.OrderedDict(
                    [
                        ("concentration", concentrations),
                        ("measurement", points),
                        ("fit", fit),
                        ("stderr", stderr),
                    ]
                )
            )
            .sort_values("concentration")
            .reset_index(drop=True)
        )


def concentrationRange(bottom, top, npoints=200, extend=0.1):
    """Logarithmically spaced concentrations for plotting.

    Useful if you want to plot a curve by fitting values to densely
    sampled points and need the concentrations at which to compute
    these points.

    Args:
        `bottom` (float)
            Lowest concentration.
        `top` (float)
            Highest concentration.
        `npoints` (int)
            Number of points.
        `extend` (float)
            After transforming to log space, extend range of points
            by this much below and above `bottom` and `top`.

    Returns:
        A numpy array of `npoints` concentrations.

    >>> numpy.allclose(concentrationRange(0.1, 100, 10, extend=0),
    ...                [0.1, 0.22, 0.46, 1, 2.15, 4.64, 10, 21.54, 46.42, 100],
    ...                atol=1e-2)
    True
    >>> numpy.allclose(concentrationRange(0.1, 100, 10),
    ...                [0.05, 0.13, 0.32, 0.79, 2.00, 5.01,
    ...                 12.59, 31.62, 79.43, 199.53],
    ...                atol=1e-2)
    True

    """
    if top <= bottom:
        raise ValueError(f"{bottom=} must be less than {top=}")
    if bottom <= 0:
        raise ValueError("{bottom=} must be greater than zero")
    if extend < 0:
        raise ValueError("{extend=} must be >= 0")

    logbottom = math.log10(bottom)
    logtop = math.log10(top)
    logrange = logtop - logbottom
    assert logrange > 0
    bottom = logbottom - logrange * extend
    top = logtop + logrange * extend

    return numpy.logspace(bottom, top, npoints)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
