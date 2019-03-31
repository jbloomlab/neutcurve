"""
==================
hillcurve
==================

Defines :class:`HillCurve` for fitting neutralization curves.
"""

import math
import collections

import pandas as pd
import scipy
import scipy.optimize


class HillCurve:
    r"""A fitted Hill curve, optionally with free baselines.

    Fits :math:`f(c) = b + \frac{t - b}{1 + (c/m)^s}`
    where :math:`f(c)` is the fraction infectivity remaining
    at concentration :math:`c`, :math:`m` is the midpoint
    of the neutralization curve, :math:`t` is the top
    value (e.g., 1), :math:`b` is the bottom value (e.g., 0),
    and :math:`s` is the slope of the curve. Because
    :math:`f(c)` is the fraction infectivity remaining, we expect
    :math:`f(c)` to get smaller as :math:`c` gets larger.
    This should lead us to fit :math:`s > 0`.

    When :math:`t = 1` and :math:`b = 0`, this equation is identical to the
    `Hill curve <https://en.wikipedia.org/wiki/Hill_equation_(biochemistry)>`_,
    except that we are calcuating the fraction **unbound** rather than
    the fraction bound.

    Args:
        `cs` (array-like)
            Concentrations of antibody / serum.
        `fs` (array-like)
            Fraction infectivity remaining at each concentration.
        `fixbottom` (bool or a float)
            If `True`, fix bottom of curve to to this value; otherwise fit.
        `fixtop` (`False` or a float)
            If `True`, fix top of curve to this value; otherwise fit.

    Attributes:
        `cs` (numpy array)
            Concentrations, sorted from low to high.
        `fs` (numpy array)
            Fraction infectivity, ordered to match sorted concentrations.
        `bottom` (float)
            Bottom of curve, :math:`b` in equation above.
        `top` (float)
            Top of curve, :math:`t` in equation above.
        `midpoint` (float)
            Midpoint of curve, :math:`m` in equation above. Note
            that the midpoint may **not** be the same as the :meth:`ic50`
            if :math:`t \ne 1` or :math:`b \ne 0`.
        `slope` (float)
            Hill slope of curve, :math:`s` in equation above.

    Use the :meth:`ic50` method to get the fitted IC50.

    Here are some examples. First, we import the modules we use, including
    `plotnine <https://plotnine.readthedocs.io>`_ for plotting:

    .. nbplot::

        >>> import scipy
        >>> import plotnine as p9
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

        >>> neut = HillCurve(cs, fs)
        >>> scipy.allclose(neut.midpoint, m)
        True
        >>> scipy.allclose(neut.slope, s)
        True
        >>> scipy.allclose(neut.top, t)
        True
        >>> scipy.allclose(neut.bottom, b)
        True

    Since we fit the curve to simulated data where the bottom was
    0.1 rather than 0, the midpoint and IC50 are different. Specifically,
    the IC50 is larger than the midpoint as you have to go past the midpoint
    to get down to 0.5 fraction infectivity:

    .. nbplot::

        >>> neut.ic50() > neut.midpoint
        True
        >>> scipy.allclose(neut.ic50(), 0.0337385586)
        True
        >>> scipy.allclose(0.5, neut.fracinfectivity(neut.ic50()))
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
        >>> neut2 = HillCurve(cs, fs2, fixbottom=b2)
        >>> scipy.allclose(neut2.midpoint, m)
        True
        >>> scipy.allclose(neut2.ic50(), m)
        True

    Now let's fit to concentrations that are all **less**
    than the midpoint, so that we never get to complete neutralization.
    The estimated IC50 is unreliable, and so will be as `None` unless
    we call :meth:`ic50` with `method` set to 'bound':

    .. nbplot::

        >>> cs3 = [1e-5 * 2**x for x in range(7)]
        >>> (cs3[-1] < m)
        True
        >>> fs3 = [HillCurve.evaluate(c, m, s, b2, t2) for c in cs3]
        >>> neut3 = HillCurve(cs3, fs3, fixbottom=b2)
        >>> neut3.ic50() is None
        True
        >>> scipy.allclose(neut3.ic50(method='bound'), cs3[-1])
        True

    We can use the :meth:`dataframe` method to get the measured
    data and fit data at selected points. First, we do this
    just at measured points:

    .. nbplot::

        >>> neut.dataframe('measured').round(3)
           concentration  measurement    fit
        0          0.002        0.995  0.995
        1          0.004        0.981  0.981
        2          0.008        0.932  0.932
        3          0.016        0.791  0.791
        4          0.032        0.522  0.522
        5          0.064        0.272  0.272
        6          0.128        0.154  0.154
        7          0.256        0.115  0.115
        8          0.512        0.104  0.104

    Then we add in one more point:

    .. nbplot::

        >>> neut.dataframe([0.6]).round(3)
           concentration  measurement    fit
        0          0.002        0.995  0.995
        1          0.004        0.981  0.981
        2          0.008        0.932  0.932
        3          0.016        0.791  0.791
        4          0.032        0.522  0.522
        5          0.064        0.272  0.272
        6          0.128        0.154  0.154
        7          0.256        0.115  0.115
        8          0.512        0.104  0.104
        9          0.600          NaN  0.103

    In reality, you'd typically just call :meth:`dataframe` with
    the default argument of 'auto' to get a good range to plot:

    .. nbplot::

        >>> df = neut.dataframe()
        >>> p = (p9.ggplot(df, p9.aes(x='concentration')) +
        ...      p9.geom_line(p9.aes(y='fit'), color=CBPALETTE[2]) +
        ...      p9.geom_point(p9.aes(y='measurement'), na_rm=True,
        ...                    color=CBPALETTE[1], size=3) +
        ...      p9.scale_x_log10(name="concentration") +
        ...      p9.scale_y_continuous(name="infectivity remaining",
        ...                            limits=(0, 1)) +
        ...      p9.theme_bw(base_size=12) +
        ...      p9.theme(figure_size=(3.5, 2.5))
        ...      )
        >>> _ = p.draw()

    """

    def __init__(self, cs, fs, *, fixbottom=False, fixtop=1):
        """See main class docstring."""
        # get data into arrays sorted by concentration
        self.cs = scipy.array(cs)
        self.fs = scipy.array(fs)
        self.fs = self.fs[self.cs.argsort()]
        self.cs = self.cs[self.cs.argsort()]

        # make initial guess for slope to have the right sign
        if self.fs[0] >= self.fs[-1]:
            self.slope = 1.5
        else:
            self.slope = -1.5

        # make initial guess for top and bottom
        if fixtop is False:
            if self.slope > 0:
                self.top = self.fs.max()
            else:
                self.top = self.fs.min()
        else:
            if not isinstance(fixtop, (int, float)):
                raise ValueError('`fixtop` is not `False` or a number')
            self.top = fixtop
        if fixbottom is False:
            if self.slope > 0:
                self.bottom = self.fs.min()
            else:
                self.bottom = self.fs.max()
        else:
            if not isinstance(fixbottom, (int, float)):
                raise ValueError('`fixbottom` is not `False` or a number')
            self.bottom = fixbottom

        # make initial guess for midpoint
        midval = (self.top - self.bottom) / 2.0
        if (self.fs > midval).all():
            if self.slope > 0:
                self.midpoint = self.cs[-1]
            else:
                self.midpoint = self.cs[0]
        elif (self.fs <= midval).all():
            if self.slope > 0:
                self.midpoint = self.cs[0]
            else:
                self.midpoint = self.cs[-1]
        else:
            # get first index where f crosses midpoint
            i = scipy.argmax((self.fs > midval)[:-1] !=
                             (self.fs > midval)[1:])
            assert (self.fs[i] > midval) != (self.fs[i + 1] > midval)
            self.midpoint = (self.cs[i] + self.cs[i + 1]) / 2.0

        # set up function and initial guesses
        if fixtop is False and fixbottom is False:
            initguess = [self.midpoint, self.slope, self.bottom, self.top]
            func = self.evaluate
        elif fixtop is False:
            initguess = [self.midpoint, self.slope, self.top]

            def func(c, m, s, t):
                return self.evaluate(c, m, s, self.bottom, t)

        elif fixbottom is False:
            initguess = [self.midpoint, self.slope, self.bottom]

            def func(c, m, s, b):
                return self.evaluate(c, m, s, b, self.top)
        else:
            initguess = [self.midpoint, self.slope]

            def func(c, m, s):
                return self.evaluate(c, m, s, self.bottom, self.top)

        (popt, pcov) = scipy.optimize.curve_fit(
                func,
                self.cs,
                self.fs,
                initguess
                )

        if fixtop is False and fixbottom is False:
            (self.midpoint, self.slope, self.top, self.bottom) = popt
        elif fixtop is False:
            (self.midpoint, self.slope, self.top) = popt
        elif fixbottom is False:
            (self.midpoint, self.slope, self.bottom) = popt
        else:
            (self.midpoint, self.slope) = popt

    def ic50(self, method='interpolate'):
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
                    return upper or lower depending on whether IC50 is
                    above or below range of concentrations. Based on
                    assumption infectivity decreases with concentration.

        Returns:
            Number giving IC50 or `None` (depending on value of `method`).

        """
        if self.top < 0.5 and self.bottom < 0.5:
            bound = 'bottom'
        elif self.top >= 0.5 and self.bottom >= 0.5:
            bound = 'upper'
        else:
            ic50 = (self.midpoint * ((self.top - 0.5) /
                    (0.5 - self.bottom))**(1.0 / self.slope))
            if (self.cs[0] <= ic50 <= self.cs[-1]):
                return ic50
            elif ic50 < self.cs[0]:
                bound = 'bottom'
            else:
                bound = 'upper'

        if method == 'bound':
            if bound == 'upper':
                return self.cs[-1]
            elif bound == 'lower':
                return self.cs[0]
            else:
                raise ValueError(f"invalid `bound` {bound}")
        elif method != 'interpolate':
            raise ValueError(f"invalid `method` of {method}")

    def fracinfectivity(self, c):
        """Fraction infectivity at `c` for fitted parameters."""
        return self.evaluate(c, self.midpoint, self.slope,
                             self.bottom, self.top)

    @staticmethod
    def evaluate(c, m, s, b, t):
        r"""Get :math:`f(c) = b + \frac{t - b}{1 + (c/m)^s}`."""
        return b + (t - b) / (1 + (c / m)**s)

    def dataframe(self, concentrations='auto'):
        """Get data frame with curve data for plotting.

        Useful if you want to get both the points and the fit
        curve to plot.

        Args:
            `concentrations` (array-like or 'auto' or 'measured')
                Concentrations for which we compute the fit values.
                If 'auto' the automatically computed from data
                range using :func:`concentrationRange`. If
                'measured' then only include measured values.

        Returns:
            A pandas DataFrame with the following columns:

              - 'concentration': concentration
              - 'fit': curve fit value at this point
              - 'measurement': value of measurement at this point,
                or numpy.nan if no measurement here.

        """
        if concentrations == 'auto':
            concentrations = concentrationRange(self.cs[0], self.cs[-1])
        elif concentrations == 'measured':
            concentrations = []
        concentrations = scipy.concatenate([self.cs, concentrations])
        n = len(concentrations)

        points = scipy.concatenate(
                [self.fs,
                 scipy.full(n - len(self.fs), scipy.nan)])

        fit = scipy.array([self.fracinfectivity(c) for c in concentrations])

        return (pd.DataFrame.from_dict(
                    collections.OrderedDict(
                        [('concentration', concentrations),
                         ('measurement', points),
                         ('fit', fit),
                         ])
                    )
                .sort_values('concentration')
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

    >>> scipy.allclose(concentrationRange(0.1, 100, 10, extend=0),
    ...                [0.1, 0.22, 0.46, 1, 2.15, 4.64, 10, 21.54, 46.42, 100],
    ...                atol=1e-2)
    True
    >>> scipy.allclose(concentrationRange(0.1, 100, 10),
    ...                [0.05, 0.13, 0.32, 0.79, 2.00, 5.01,
    ...                 12.59, 31.62, 79.43, 199.53],
    ...                atol=1e-2)
    True

    """
    if top <= bottom:
        raise ValueError('`bottom` must be less than `top`')
    if bottom <= 0:
        raise ValueError('`bottom` must be greater than zero')
    if extend < 0:
        raise ValueError('`extend` must be >= 0')

    logbottom = math.log10(bottom)
    logtop = math.log10(top)
    logrange = logtop - logbottom
    assert logrange > 0
    bottom = logbottom - logrange * extend
    top = logtop + logrange * extend

    return scipy.logspace(bottom, top, npoints)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
