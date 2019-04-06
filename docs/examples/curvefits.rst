.. _curvefits_example:

Fitting curves to real data
===========================

.. contents:: Contents
   :local:

.. code-links:: full
    :timeout: 200

Overview
---------

When you have a real dataset, the easiest way to fit it is via the
:class:`neutcurve.curvefits.CurveFits` class, which
streamlines the fitting and plotting of a set of
:class:`neutcurve.hillcurve.HillCurve` neutralization curves.
Here we demonstrate how to do this.

Importing the necessary packages
--------------------------------

First, import the necessary Python modules. In addition to the
``neutcurve`` package itself, we also use
`pandas <https://pandas.pydata.org/>`__ to hold the data:

.. nbplot::

    >>> import urllib.request
    >>>
    >>> import pandas as pd
    >>>
    >>> import neutcurve
    >>> from neutcurve.colorschemes import CBMARKERS, CBPALETTE

Set some pandas display options:

.. nbplot::

    >>> pd.set_option('display.float_format', '{:.3g}'.format)
    >>> pd.set_option('display.max_columns', 20)
    >>> pd.set_option('display.width', 500)

Preparing the data
------------------

Now we get example data to fit. We use as our example the neutralization
of variants A/WSN/1933 (H1N1) influenza by the broadly neutralizing
antibody FI6v3 and strain-specific antibody H17-L19 from `Fig 6a,b of
Doud et al
(2018) <https://www.nature.com/articles/s41467-018-03665-3#Fig6>`__.
These data in numerical form in a CSV file are `available
here <https://github.com/jbloomlab/neutcurve/blob/master/docs/examples/data/Doud_et_al_2018-neutdata.csv>`__.
We start by downloading these data and reading them into a pandas
DataFrame:

.. nbplot::

    >>> fi6v3_datafile = ('https://raw.githubusercontent.com/jbloomlab/neutcurve/'
    ...                   'master/docs/examples/data/Doud_et_al_2018-neutdata.csv')
    ...
    >>> with urllib.request.urlopen(fi6v3_datafile) as f:
    ...     data = pd.read_csv(f)

Here are the first few lines of the data frame:

.. nbplot::

    >>> data.head()
       serum virus  replicate  concentration  fraction infectivity
    0  FI6v3    WT          1       0.000205                  1.01
    1  FI6v3    WT          1       0.000478                 0.942
    2  FI6v3    WT          1        0.00112                 0.993
    3  FI6v3    WT          1         0.0026                 0.966
    4  FI6v3    WT          1        0.00607                 0.957

And here are the last few lines:

.. nbplot::

    >>> data.tail()
           serum  virus  replicate  concentration  fraction infectivity
    427  H17-L19  V135T          3          0.386                  1.02
    428  H17-L19  V135T          3            0.9                     1
    429  H17-L19  V135T          3            2.1                 0.959
    430  H17-L19  V135T          3            4.9                 0.991
    431  H17-L19  V135T          3           11.4                 0.747

As can be seen above, the data are organized into five columns, all of
which must be present. These columns are: 
  - *serum*: the name of the
    serum (or antibody). FI6v3 and H17-L19 are actually antibodies, not
    sera, but :class:`neutcurve.curvefits.CurveFits` is set up to refer to
    things as serum. 
  - *virus*: the name of the virus being neutralized by
    the serum. 
  - *replicate*: the replicate label for the measurement.
    Although you can have just one replicate, it’s good experimental
    practice to have several. All the replicates for a given virus / serum
    combination must have been measured at the same concentrations. 
  - *concentration*: the concentration of the serum. 
  - *fraction infectivity*: the fraction infectivity of the virus at this
    concentration of the serum measured in this replicate.

Note that the data are in `tidy form <https://cran.r-project.org/web/packages/tidyr/vignettes/tidy-data.html>`__;
you must make your data frame tidy before you can analyze it with
:class:`neutcurve.curvefits.CurveFits`.

Fitting the curves
------------------

Once you have the tidy data frame, it’s easy to pass it to
:class:`neutcurve.curvefits.CurveFits`. We expect all of these
antibodies to go to complete neutralization when they are effective, so
we use the `fixbottom=0` argument (see
:class:`neutcurve.hillcurve.HillCurves` and :ref:`hillcurve_example` for more details about the
`fixtop` and `fixbottom` options):

.. nbplot::

    >>> fits = neutcurve.CurveFits(data,
    ...                            fixbottom=0,
    ...                            )

Now we can look at the different sera for which we have fit curves:

.. nbplot::

    >>> fits.sera
    ['FI6v3', 'H17-L19']

We can also look at the viruses measured against each serum:

.. nbplot::

    >>> for serum in fits.sera:
    ...     print(f"Viruses measured against {serum}:\n" +
    ...           str(fits.viruses[serum]))
    Viruses measured against FI6v3:
    ['WT', 'K(-8T)', 'P80D', 'V135T', 'K280A', 'K280S', 'K280T', 'N291S', 'M17L-HA2', 'G47R-HA2']
    Viruses measured against H17-L19:
    ['WT', 'V135T']

We can also look at the replicates for each serum and virus. Here we
just do that for serum *FI6v3* and virus *WT*. See how in addition to
the three replicates we have passed, there is also now an “average”
replicate that is automatically computed from the average of the other
replicates:

.. nbplot::

    >>> fits.replicates[('FI6v3', 'WT')]
    ['1', '2', '3', 'average']

Looking at a specific curve
---------------------------

We can use the :meth:`neutcurve.curvefits.CurveFits.getCurve` method
to get the :class:`neutcurve.hillcurve.HillCurve` that was fit for a
particular serum / virus / replicate combination. For instance, here we
do that for *serum* FI6v3 versus *virus* WT for replicate *1*. We then
plot the curve and print the IC50:

.. nbplot::

    >>> curve = fits.getCurve(serum='FI6v3', virus='WT', replicate='1')
    >>> print(f"The IC50 is {curve.ic50():.3g}")
    The IC50 is 0.0167
    >>> fig, ax = curve.plot()

:class:`neutcurve.curvefits.CurveFits` also calculates the average and
standard error of the measurements for each serum / virus, and fits them
under a replicate name of “average”. Here is the fit to the average of
the data for *serum* FI6v3 and *virus* WT. Note how the plot now also
shows error bars indicating the standard error:

.. nbplot::

    >>> curve = fits.getCurve(serum='FI6v3', virus='WT', replicate='average')
    >>> print(f"The IC50 is {curve.ic50():.3g}")
    The IC50 is 0.0195
    >>> fig, ax = curve.plot()

Accessing the fit parameters for all curves
-------------------------------------------

You can get the fit parameters for the curves using
:meth:`neutcurve.curvefits.CurveFits.fitParams`. By default, this just
gets the fits for the average of the replicates. The parameters are all
of those fit by a :class:`neutcurve.hillcurve.HillCurve`, plus the
IC50 in several forms to accurately represent interpolated IC50s (IC50
within range of data) versus IC50s where we can just get the bound from
the upper or lower limits of the data:

.. nbplot::

    >>> fits.fitParams()
          serum     virus replicate  nreplicates   ic50    ic50_bound ic50_str  midpoint  slope  top  bottom
    0     FI6v3        WT   average            3 0.0195  interpolated   0.0195    0.0195      3    1       0
    1     FI6v3    K(-8T)   average            3 0.0325  interpolated   0.0325    0.0325    2.9    1       0
    2     FI6v3      P80D   average            3 0.0124  interpolated   0.0124    0.0124    2.3    1       0
    3     FI6v3     V135T   average            3 0.0241  interpolated   0.0241    0.0241   2.09    1       0
    4     FI6v3     K280A   average            3 0.0142  interpolated   0.0142    0.0142   3.18    1       0
    5     FI6v3     K280S   average            3 0.0389  interpolated   0.0389    0.0389    2.9    1       0
    6     FI6v3     K280T   average            3 0.0392  interpolated   0.0392    0.0392   2.33    1       0
    7     FI6v3     N291S   average            3  0.106  interpolated    0.106     0.106   2.77    1       0
    8     FI6v3  M17L-HA2   average            3  0.022  interpolated    0.022     0.022   2.69    1       0
    9     FI6v3  G47R-HA2   average            3 0.0356  interpolated   0.0356    0.0356   3.32    1       0
    10  H17-L19        WT   average            3   0.11  interpolated     0.11      0.11   4.76    1       0
    11  H17-L19     V135T   average            3   11.4         upper    >11.4      15.5   2.77    1       0

Looking above, you can see how the IC50 is handled depending on if it is
interpolated (in the range of concentrations used in the experiments)
versus outside the range of concentrations. In the table above, all of
the IC50s are interpolated **except** the last row (H17-L19 versus
V135T), which is just provided as an upper bound equal to the highest
concentration used in the experiment (the actual IC50 is greater than
this upper bound). We do **not** attempt to extrapolate IC50s outside
the data range as this is unreliable.

Note that by default, :meth:`neutcurve.curvefits.CurveFits.fitParams`
only returns the fitted params for the averages, as in the above table.
If you want to also return them for individual replicates, using the
`average_only=False` argument. Here we do this, showing only the first
few entries in the returned data frame; now there are now values for
each replicate as well as the average of replicates:

.. nbplot::

    >>> fits.fitParams(average_only=False).head()
       serum   virus replicate  nreplicates   ic50    ic50_bound ic50_str  midpoint  slope  top  bottom
    0  FI6v3      WT         1          NaN 0.0167  interpolated   0.0167    0.0167    2.5    1       0
    1  FI6v3      WT         2          NaN  0.019  interpolated    0.019     0.019   2.51    1       0
    2  FI6v3      WT         3          NaN 0.0152  interpolated   0.0152    0.0152   1.88    1       0
    3  FI6v3      WT   average            3 0.0195  interpolated   0.0195    0.0195      3    1       0
    4  FI6v3  K(-8T)         1          NaN 0.0308  interpolated   0.0308    0.0308   2.62    1       0

Note that the “average” is the curve fit to the average of the data
points, not the average of the fit parameters for individual curves.

Plotting the curves
-------------------

One of the most useful feature of
:class:`neutcurve.curvefits.CurveFits` are that they have methods to
easily generate multi-panel plots of the curves.

Plotting each serum against multiple viruses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Often you will have measured each serum against several different viral
variants. You can then plot these curves using
:meth:`neutcurve.curvefits.CurveFits.plotSera` as below:

.. nbplot::

    >>> fig, axes = fits.plotSera(xlabel='concentration (ug/ml)')

The above plot attempts to put all the viruses measured against each
serum on the same subplot, but is cognizant of the fact that it becomes
uninterpretable if there are too many viruses on the same plot.
Therefore, it only shows a maximum of `max_viruses_per_subplot` (which
by default is 5) curves per subplot.

In fact, that is still perhaps too many curves per plot for this data set. So we can
customize the plot by adjusting that parameter. Below we adjust to just
four viruses per subplot, and also use `ncol=2` to specify that we
want two columns:

.. nbplot::

    >>> fig, axes = fits.plotSera(max_viruses_per_subplot=4,
    ...                           ncol=2,
    ...                           xlabel='concentration (ug/ml)')




The above plots all have a different legend for each subplot. This is
necessary because the number of different viruses being plotted exceeds
the numbers of colors / markers specified to
:meth:`neutcurve.curvefits.CurveFits.plotSera` via its `colors` and
`markers` arguments, so there aren’t enough colors / markers to give
each virus a unique one.

However, if we reduce the number of viruses we are showing, we then get
a nice shared legend. Here we do this, using the `viruses` argument to
specify that we just show some of the viruses:

.. nbplot::

    >>> fig, axes = fits.plotSera(viruses=['WT', 'N291S', 'K280S', 'V135T'],
    ...                           xlabel='concentration (ug/ml)') 




Similar to how the above plot uses the `viruses` argument to plot just
some viruses, we can also use the `sera` argument to plot just some of
the sera (in this case, just H17-L19):

.. nbplot::

    >>> fig, axes = fits.plotSera(sera=['H17-L19'],
    ...                           xlabel='concentration (ug/ml)')




There are various additional options to
:meth:`neutcurve.curvefits.CurveFits.plotSera` that can further
fine-tune the plots; see the docstring for that method for more details.

Plotting each replicate
~~~~~~~~~~~~~~~~~~~~~~~

Another type of plot that is sometimes useful is one that shows all the
replicates for each serum / virus combination. Such a plot is easily
generated using :meth:`neutcurve.curvefits.CurveFits.plotReplicates`
as below:

.. nbplot::

    >>> fig, axes = fits.plotReplicates(xlabel='concentration (ug/ml)',
    ...                                 legendtitle='replicate')




See the method docstring for
:meth:`neutcurve.curvefits.CurveFits.plotReplicates` for ways to
further customize these plots.

Plotting average for each serum-virus combination
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Another type of plot that is useful is one that simply shows the
replicate-average for each serum-virus on its own subplot. This plot can
be generated by again using
:meth:`neutcurve.curvefits.CurveFits.plotReplicates` but now
specifying that we only show the average value by setting
`average_only=True`. Below we do that, also using `colors` to
specify that the single curve for each subplot is black:

.. nbplot::

    >>> fig, axes = fits.plotReplicates(xlabel='concentration (ug/ml)',
    ...                                 average_only=True,
    ...                                 colors=['black'])




Customized plot arrangements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are obviously many other ways that it’s possible to lay out the
different curves for sera / viruses / replicates on subplots. You can
make an arbitrarily customized layout using
:meth:`neutcurve.curvefits.CurveFits.plotGrid` where you explicitly
pass the curves to put at each subplot in the plot.

Below we illustrate how to do this to create a plot that essentially
mimics what is shown in `Fig 6a,b of Doud et al
(2018) <https://www.nature.com/articles/s41467-018-03665-3#Fig6>`__
(although those published plots were not generated using this program).
Note that in doing this below, we use the colors and markers defined by
`CBPALETTE` and `CBMARKERS` in :mod:`neutcurve.colorschemes`:

.. nbplot::

    >>> fig, axes = fits.plotGrid(
    ...                 {
    ...                  # upper right: FI6v3 versus WT, K280S, K280T, K280A
    ...                  (0, 0): ('FI6v3',
    ...                           [{'serum': 'FI6v3', 'virus': 'WT',
    ...                             'replicate': 'average', 'color': CBPALETTE[0],
    ...                             'marker': CBMARKERS[0], 'label':'WT'},
    ...                            {'serum': 'FI6v3', 'virus': 'K280S',
    ...                             'replicate': 'average', 'color': CBPALETTE[1],
    ...                             'marker': CBMARKERS[1], 'label':'K280S'},
    ...                            {'serum': 'FI6v3', 'virus': 'K280T',
    ...                             'replicate': 'average', 'color': CBPALETTE[2],
    ...                             'marker': CBMARKERS[2], 'label':'K280T'},
    ...                            {'serum': 'FI6v3', 'virus': 'K280A',
    ...                             'replicate': 'average', 'color': CBPALETTE[3],
    ...                             'marker': CBMARKERS[3], 'label':'K280A'},
    ...                            ]
    ...                           ),
    ...                  # upper center: FI6v3 versus WT, N291S
    ...                  (0, 1): ('FI6v3',
    ...                           [{'serum': 'FI6v3', 'virus': 'WT',
    ...                             'replicate': 'average', 'color': CBPALETTE[0],
    ...                             'marker': CBMARKERS[0], 'label': 'WT'},
    ...                            {'serum': 'FI6v3', 'virus': 'N291S',
    ...                             'replicate': 'average', 'color': CBPALETTE[1],
    ...                             'marker': CBMARKERS[1], 'label': 'N291S'},
    ...                            ]
    ...                           ),
    ...                  # upper right: FI6v3 versus WT, G47R-HA2
    ...                  (0, 2): ('FI6v3',
    ...                           [{'serum': 'FI6v3', 'virus': 'WT',
    ...                             'replicate': 'average', 'color': CBPALETTE[0],
    ...                             'marker': CBMARKERS[0], 'label': 'WT'},
    ...                            {'serum': 'FI6v3', 'virus': 'G47R-HA2',
    ...                             'replicate': 'average', 'color': CBPALETTE[1],
    ...                             'marker': CBMARKERS[1], 'label': 'G47R(HA2)'},
    ...                            ]
    ...                           ),
    ...                  # middle right: FI6v3 versus WT, K(-8T)
    ...                  (1, 0): ('FI6v3',
    ...                           [{'serum': 'FI6v3', 'virus': 'WT',
    ...                             'replicate': 'average', 'color': CBPALETTE[0],
    ...                             'marker': CBMARKERS[0], 'label': 'WT'},
    ...                            {'serum': 'FI6v3', 'virus': 'K(-8T)',
    ...                             'replicate': 'average', 'color': CBPALETTE[1],
    ...                             'marker': CBMARKERS[1], 'label': 'K(-8T)'},
    ...                            ]
    ...                           ),
    ...                  # middle center: FI6v3 versus WT, M17L-HA2
    ...                  (1, 1): ('FI6v3',
    ...                           [{'serum': 'FI6v3', 'virus': 'WT',
    ...                             'replicate': 'average', 'color': CBPALETTE[0],
    ...                             'marker': CBMARKERS[0], 'label': 'WT'},
    ...                            {'serum': 'FI6v3', 'virus': 'M17L-HA2',
    ...                             'replicate': 'average', 'color': CBPALETTE[1],
    ...                             'marker': CBMARKERS[1], 'label': 'M17L(HA2)'},
    ...                            ]
    ...                           ),
    ...                  # middle right: FI6v3 versus WT, P80D, V135T
    ...                  (1, 2): ('FI6v3',
    ...                           [{'serum': 'FI6v3', 'virus': 'WT',
    ...                             'replicate': 'average', 'color': CBPALETTE[0],
    ...                             'marker': CBMARKERS[0], 'label': 'WT'},
    ...                            {'serum': 'FI6v3', 'virus': 'P80D',
    ...                             'replicate': 'average', 'color': CBPALETTE[1],
    ...                             'marker': CBMARKERS[1], 'label': 'P80D'},
    ...                            {'serum': 'FI6v3', 'virus': 'V135T',
    ...                             'replicate': 'average', 'color': CBPALETTE[2],
    ...                             'marker': CBMARKERS[2], 'label': 'V135T'},
    ...                            ]
    ...                           ),
    ...                  # middle left: H17-L19 versus WT, V135T
    ...                  (2, 0): ('H17-L19',
    ...                           [{'serum': 'H17-L19', 'virus': 'WT',
    ...                             'replicate': 'average', 'color': CBPALETTE[0],
    ...                             'marker': CBMARKERS[0], 'label': 'WT'},
    ...                            {'serum': 'H17-L19', 'virus': 'V135T',
    ...                             'replicate': 'average', 'color': CBPALETTE[2],
    ...                             'marker': CBMARKERS[1], 'label': 'V135T'},
    ...                            ]
    ...                           ),
    ...                  },
    ...                 xlabel='concentration (ug/ml)',
    ...                 )




