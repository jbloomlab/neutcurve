"""
==================
curvefits
==================
Defines :class:`CurveFits` to fit curves and display / plot results.
"""

import collections
import itertools
import math

import dmslogo.facet
import dmslogo.utils

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import pandas as pd

import neutcurve
from neutcurve.colorschemes import CBMARKERS, CBPALETTE


class CurveFits:
    """Fit and display :class:`neutcurve.hillcurve.HillCurve` curves.

    Args:
        `data` (pandas DataFrame)
            Tidy dataframe with data.
        `conc_col` (str)
            Column in `data` with concentrations of serum.
        `fracinf_col` (str)
            Column in `data` with fraction infectivity.
        `serum_col` (str)
            Column in `data` with serum name.
        `virus_col` (str)
            Column in `data` with name of virus being neutralized.
        `replicate_col` (str`)
            Column in data with name of replicate of this measurement.
            Replicates must all have the same concentrations for each
            serum / virus combination. Replicates can **not** be named
            'average' as we compute the average from the replicates.
        `fixbottom` (`False` or float)
            Same meaning as for :class:`neutcurve.hillcurve.HillCurve`.
        `fixtop` (`False` or float)
            Same meaning as for :class:`neutcurve.hillcurve.HillCurve`.
        `infectivity_or_neutralized` ({'infectivity', 'neutralized'})
            Same meaning as for :class:`neutcurve.hillcurve.HillCurve`.

    Attributes of a :class:`CurveFits` include all args except `data` plus:
        `df` (pandas DataFrame)
            Copy of `data` that only has relevant columns, has additional rows
            with `replicate_col` of 'average' that hold replicate averages, and
            added columns 'stderr' (standard error of fraction infectivity
            for 'average' if multiple replicates, otherwise `nan`).
        `sera` (list)
            List of all serum names in `serum_col` of `data`, in order
            they occur in `data`.
        `viruses` (dict)
            For each serum in `sera`, `viruses[serum]` gives all viruses
            for that serum in the order they occur in `data`.
        `replicates` (dict)
            `replicates[(serum, virus)]` is list of all replicates for
            that serum and virus in the order they occur in `data`.
        `allviruses` (list)
            List of all viruses.

    """

    # names commonly used for wildtype virus
    _WILDTYPE_NAMES = ('WT', 'wt', 'wildtype', 'Wildtype',
                       'wild type', 'Wild type')

    def __init__(self,
                 data,
                 *,
                 conc_col='concentration',
                 fracinf_col='fraction infectivity',
                 serum_col='serum',
                 virus_col='virus',
                 replicate_col='replicate',
                 infectivity_or_neutralized='infectivity',
                 fixbottom=0,
                 fixtop=1,
                 ):
        """See main class docstring."""
        # make args into attributes
        self.conc_col = conc_col
        self.fracinf_col = fracinf_col
        self.serum_col = serum_col
        self.virus_col = virus_col
        self.replicate_col = replicate_col
        self.fixbottom = fixbottom
        self.fixtop = fixtop
        self._infectivity_or_neutralized = infectivity_or_neutralized

        # check for required columns
        cols = [self.serum_col, self.virus_col, self.replicate_col,
                self.conc_col, self.fracinf_col]
        if len(cols) != len(set(cols)):
            raise ValueError('duplicate column names:\n\t' + '\n\t'.join(cols))
        if not (set(cols) <= set(data.columns)):
            raise ValueError('`data` lacks required columns, which are:\n\t' +
                             '\n\t'.join(cols))

        # create `self.df`, ensure that replicates are str rather than number
        self.df = (data[cols]
                   .assign(**{replicate_col: lambda x: (x[replicate_col]
                                                        .astype(str))
                              })
                   )

        # create sera / viruses / replicates attributes, error check them
        self.sera = self.df[self.serum_col].unique().tolist()
        self.viruses = {}
        self.replicates = {}
        for serum in self.sera:
            serum_data = self.df.query(f"{self.serum_col} == @serum")
            serum_viruses = serum_data[self.virus_col].unique().tolist()
            self.viruses[serum] = serum_viruses
            for virus in serum_viruses:
                virus_data = serum_data.query(f"{self.virus_col} == @virus")
                virus_reps = virus_data[self.replicate_col].unique().tolist()
                if 'average' in virus_reps:
                    raise ValueError('A replicate is named "average". This is '
                                     'not allowed as that name is used for '
                                     'replicate averages.')
                self.replicates[(serum, virus)] = virus_reps + ['average']
                for i, rep1 in enumerate(virus_reps):
                    conc1 = (virus_data
                             .query(f"{self.replicate_col} == @rep1")
                             [self.conc_col]
                             .sort_values()
                             .tolist()
                             )
                    if len(conc1) != len(set(conc1)):
                        raise ValueError('duplicate concentrations for '
                                         f"{serum}, {virus}, {rep1}")
                    for rep2 in virus_reps[i + 1:]:
                        conc2 = (virus_data
                                 .query(f"{self.replicate_col} == @rep1")
                                 [self.conc_col]
                                 .sort_values()
                                 .tolist()
                                 )
                        if conc1 != conc2:
                            raise ValueError(f"replicates {rep1} and {rep2} "
                                             'have different concentrations '
                                             f"for {serum}, {virus}")
        self.allviruses = collections.OrderedDict()
        for serum in self.sera:
            for virus in self.viruses[serum]:
                self.allviruses[virus] = True
        self.allviruses = list(self.allviruses.keys())

        if pd.isnull(self.allviruses).any():
            raise ValueError(f"a virus has name NaN:\n{self.allviruses}")
        if pd.isnull(self.sera).any():
            raise ValueError(f"a serum has name NaN:\n{self.sera}")

        # compute replicate average and add 'stderr'
        if 'stderr' in self.df.columns:
            raise ValueError('`data` has column "stderr"')
        avg_df = (self.df
                  .groupby([self.serum_col, self.virus_col, self.conc_col],
                           observed=True)
                  [self.fracinf_col]
                  # sem is sample stderr, evaluates to NaN when just 1 rep
                  .aggregate(['mean', 'sem', 'count'])
                  .rename(columns={'mean': self.fracinf_col,
                                   'sem': 'stderr',
                                   })
                  .reset_index()
                  .assign(**{replicate_col: 'average'})
                  )
        self.df = pd.concat([self.df, avg_df],
                            ignore_index=True,
                            sort=False,
                            )

        self._hillcurves = {}  # curves computed by `getCurve` cached here
        self._fitparams = {}  # cache data frame computed by `fitParams`

    def getCurve(self, *, serum, virus, replicate):
        """Get the fitted curve for this sample.

        Args:
            `serum` (str)
                Name of a valid serum.
            `virus` (str)
                Name of a valid virus for `serum`.
            `replicate` (str)
                Name of a valid replicate for `serum` and `virus`, or
                'average' for the average of all replicates.

        Returns:
            A :class:`neutcurve.hillcurve.HillCurve`.

        """
        key = (serum, virus, replicate)

        if key not in self._hillcurves:
            if serum not in self.sera:
                raise ValueError(f"invalid `serum` of {serum}")
            if virus not in self.viruses[serum]:
                raise ValueError(f"invalid `virus` of {virus} for "
                                 f"`serum` of {serum}")
            if replicate not in self.replicates[(serum, virus)]:
                raise ValueError(f"invalid `replicate` of {replicate} for "
                                 f"`serum` of {serum} and `virus` of {virus}")

            idata = self.df.query(f"({self.serum_col} == @serum) & "
                                  f"({self.virus_col} == @virus) & "
                                  f"({self.replicate_col} == @replicate)")
            if len(idata) < 1:
                raise RuntimeError(f"no data for serum {serum} virus {virus}")

            if idata['stderr'].isna().all():
                fs_stderr = None
            elif idata['stderr'].isna().any():
                raise RuntimeError('`stderr` has only some entries NaN\n' +
                                   str(idata))
            else:
                fs_stderr = idata['stderr']

            try:
                curve = neutcurve.HillCurve(
                        cs=idata[self.conc_col],
                        fs=idata[self.fracinf_col],
                        fs_stderr=fs_stderr,
                        fixbottom=self.fixbottom,
                        fixtop=self.fixtop,
                        infectivity_or_neutralized=(
                            self._infectivity_or_neutralized)
                        )
            except RuntimeError as e:
                idata.to_csv('temp.csv', index=False)
                # following here: https://stackoverflow.com/a/46091127
                raise RuntimeError('Error while fitting HillCurve for '
                                   f"serum {serum}, virus {virus}, "
                                   f"replicate {replicate}.\nData are:\n" +
                                   str(idata)) from e

            self._hillcurves[key] = curve

        return self._hillcurves[key]

    def fitParams(self,
                  *,
                  average_only=True,
                  ics=(50,),
                  ics_precision=0,
                  ):
        """Get data frame with curve fitting parameters.

        Args:
            `average_only` (bool)
                If `True`, only get parameters for average across replicates.
            `ics` (iterable)
                Include ICXX for each number in this list, where the number
                is the percent neutralized. So if `ics` only contains 50,
                we include the IC50. If it includes 95, we include the IC95.
            `ics_precision` (int)
                Include this many digits after decimal when creating the
                ICXX columns.

        Returns:
            A pandas DataFrame with fit parameters for each serum / virus /
            replicate as defined for a :mod:`neutcurve.hillcurve.HillCurve`.
            Columns:

              - 'serum'
              - 'virus'
              - 'replicate'
              - 'nreplicates': number of replicates for average, NaN otherwise.
              - 'icXX': ICXX or its bound as a number, where `XX` is each
                number in `ics`.
              - 'icXX_bound': string indicating if ICXX interpolated from data,
                or is an upper or lower bound.
              - 'icXX_str': ICXX represented as string, with > or < indicating
                if it is an upper or lower bound.
              - 'midpoint': midpoint of curve, same as IC50 only if bottom
                and top are 0 and 1.
              - 'slope': Hill slope of curve.
              - 'top': top of curve.
              - 'bottom': bottom of curve.

        """
        ics = tuple(ics)
        ic_colprefixes = [f"ic{{:.{ics_precision}f}}".format(ic) for ic in ics]
        if len(ic_colprefixes) != len(set(ic_colprefixes)):
            raise ValueError('column names for ICXX not unique.\n'
                             'Either you have duplicate entries in `ics` '
                             'or you need to increase `ics_precision`.')

        key = (average_only, ics, ics_precision)

        if key not in self._fitparams:
            d = collections.defaultdict(list)
            params = ['midpoint', 'slope', 'top', 'bottom']
            for serum in self.sera:
                for virus in self.viruses[serum]:
                    replicates = self.replicates[(serum, virus)]
                    nreplicates = sum(r != 'average' for r in replicates)
                    assert nreplicates == len(replicates) - 1
                    if average_only:
                        replicates = ['average']
                    for replicate in replicates:
                        curve = self.getCurve(serum=serum,
                                              virus=virus,
                                              replicate=replicate
                                              )
                        d['serum'].append(serum)
                        d['virus'].append(virus)
                        d['replicate'].append(replicate)
                        if replicate == 'average':
                            d['nreplicates'].append(nreplicates)
                        else:
                            d['nreplicates'].append(float('nan'))
                        for ic, colprefix in zip(ics, ic_colprefixes):
                            f = ic / 100
                            d[colprefix].append(curve.icXX(f, method='bound'))
                            d[f"{colprefix}_bound"].append(curve.icXX_bound(f))
                            d[f"{colprefix}_str"].append(curve.icXX_str(f))
                        for param in params:
                            d[param].append(getattr(curve, param))

            ic_cols = []
            for prefix in ic_colprefixes:
                ic_cols += [prefix, f"{prefix}_bound", f"{prefix}_str"]
            self._fitparams[key] = (
                    pd.DataFrame(d)
                    [['serum', 'virus', 'replicate', 'nreplicates']
                     + ic_cols + params]
                    .assign(nreplicates=lambda x: (x['nreplicates']
                                                   .astype('Int64'))
                            )
                    )

        return self._fitparams[key]

    def plotSera(self,
                 *,
                 ncol=4,
                 nrow=None,
                 sera='all',
                 viruses='all',
                 ignore_serum_virus=None,
                 colors=CBPALETTE,
                 markers=CBMARKERS,
                 virus_to_color_marker=None,
                 max_viruses_per_subplot=5,
                 multi_serum_subplots=True,
                 all_subplots=_WILDTYPE_NAMES,
                 titles=None,
                 vlines=None,
                 **kwargs,
                 ):
        """Plot grid with replicate-average of viruses for each serum.

        Args:
            `ncol`, `nrow` (int or `None`)
                Specify one of these to set number of columns or rows,
                other should be `None`.
            `sera` ('all' or list)
                Sera to include on plot, in this order.
            `viruses` ('all' or list)
                Viruses to include on plot, in this order unless one
                is specified in `all_subplots`.
            `ignore_serum_virus` (`None` or dict)
                Specific serum / virus combinations to ignore (not plot). Key
                by serum, and then list viruses to ignore.
            `colors` (iterable)
                List of colors for different viruses.
            `markers` (iterable)
                List of markers for different viruses.
            `virus_to_color_marker` (dict or `None`)
                Optionally specify a specific color and for each virus as
                2-tuples `(color, marker)`. If you use this option, `colors`
                and `markers` are ignored.
            `max_viruses_per_subplot` (int)
                Maximum number of viruses to show on any subplot.
            `multi_serum_subplots` (bool)
                If a serum has more than `max_virus_per_subplot` viruses,
                do we make multiple subplots for it or raise an error?
            `all_subplots` (iterable)
                If making multiple subplots for serum, which viruses
                do we show on all subplots? These are also shown first.
            `titles` (`None` or list)
                Specify custom titles for each subplot different than `sera`.
            `vlines` (`None` or dict)
                Add vertical lines to plots. Keyed by serum name, values
                are lists of dicts with a key 'x' giving x-location of vertical
                line, and optional keys 'linewidth', 'color', and 'linestyle'.
            `**kwargs`
                Other keyword arguments that can be passed to
                :meth:`CurveFits.plotGrid`.

        Returns:
            The 2-tuple `(fig, axes)` of matplotlib figure and 2D axes array.

        """
        sera, viruses = self._sera_viruses_lists(sera, viruses)
        viruses = list(collections.OrderedDict.fromkeys(viruses))

        if titles is None:
            titles = sera
        elif len(sera) != len(titles):
            raise ValueError(f"`titles`, `sera` != length:\n{titles}\n{sera}")

        if max_viruses_per_subplot < 1:
            raise ValueError('`max_viruses_per_subplot` must be at least 1')

        # get color scheme for viruses
        if virus_to_color_marker:
            extra_viruses = set(viruses) - set(virus_to_color_marker.keys())
            if extra_viruses:
                raise ValueError('viruses not in `virus_to_color_marker`: ' +
                                 str(extra_viruses))
        elif len(viruses) <= min(len(colors), len(markers)):
            # can share scheme among subplots
            ordered_viruses = ([v for v in viruses if v in all_subplots] +
                               [v for v in viruses if v not in all_subplots])
            virus_to_color_marker = {v: (c, m) for (v, c, m) in
                                     zip(ordered_viruses, colors, markers)}
        elif min(len(colors), len(markers)) < max_viruses_per_subplot:
            raise ValueError('`max_viruses_per_subplot` larger than '
                             'number of colors or markers')
        else:
            virus_to_color_marker = None

        # Build a list of plots appropriate for `plotGrid`.
        # Code is complicated because we could have several curve
        # per serum, and in that case need to share viruses in
        # `all_subplots` among curves.
        plotlist = []
        vlines_list = []
        for serum, title in zip(sera, titles):
            if ignore_serum_virus and serum in ignore_serum_virus:
                ignore_virus = ignore_serum_virus[serum]
            else:
                ignore_virus = {}
            curvelist = []
            ivirus = 0
            serum_shared_viruses = [v for v in self.viruses[serum] if
                                    (v in viruses) and (v in all_subplots) and
                                    (v not in ignore_virus)]
            serum_unshared_viruses = [v for v in self.viruses[serum] if
                                      (v in viruses) and
                                      (v not in all_subplots) and
                                      (v not in ignore_virus)]
            unshared = int(bool(len(serum_unshared_viruses)))
            if len(serum_shared_viruses) > max_viruses_per_subplot - unshared:
                raise ValueError(f"serum {serum} has too many subplot-shared "
                                 'viruses (in `all_subplots`) relative to '
                                 'value of `max_viruses_per_subplot`:\n'
                                 f"{serum_shared_viruses} is more than "
                                 f"{max_viruses_per_subplot} viruses.")
            shared_curvelist = []
            for virus in serum_shared_viruses + serum_unshared_viruses:
                if ivirus >= max_viruses_per_subplot:
                    if multi_serum_subplots:
                        plotlist.append((title, curvelist))
                        if vlines and (serum in vlines):
                            vlines_list.append(vlines[serum])
                        else:
                            vlines_list.append(None)
                        curvelist = list(shared_curvelist)
                        ivirus = len(curvelist)
                        assert ivirus < max_viruses_per_subplot
                    else:
                        raise ValueError(f"serum {serum} has more than "
                                         '`max_viruses_per_subplot` viruses '
                                         'and `multi_serum_subplots` is False')
                if virus_to_color_marker:
                    color, marker = virus_to_color_marker[virus]
                else:
                    color = colors[ivirus]
                    marker = markers[ivirus]
                curvelist.append({'serum': serum,
                                  'virus': virus,
                                  'replicate': 'average',
                                  'label': virus,
                                  'color': color,
                                  'marker': marker,
                                  })
                if virus in serum_shared_viruses:
                    shared_curvelist.append(curvelist[-1])
                ivirus += 1
            if curvelist:
                plotlist.append((title, curvelist))
                if vlines and (serum in vlines):
                    vlines_list.append(vlines[serum])
                else:
                    vlines_list.append(None)
        if not plotlist:
            raise ValueError('no curves for these sera / viruses')

        # get number of columns
        if (nrow is not None) and (ncol is not None):
            raise ValueError('either `ncol` or `nrow` must be `None`')
        elif isinstance(nrow, int) and nrow > 0:
            ncol = math.ceil(len(plotlist) / nrow)
        elif not (isinstance(ncol, int) and ncol > 0):
            raise ValueError('`nrow` or `ncol` must be integer > 0')

        # convert plotlist to plots dict for `plotGrid`
        plots = {}
        vlines_axkey = {}
        assert len(plotlist) == len(vlines_list)
        for iplot, (plot, ivline) in enumerate(zip(plotlist, vlines_list)):
            irow = iplot // ncol
            icol = iplot % ncol
            plots[(irow, icol)] = plot
            if ivline:
                vlines_axkey[(irow, icol)] = ivline

        if virus_to_color_marker and 'orderlegend' not in kwargs:
            orderlegend = virus_to_color_marker.keys()
        else:
            orderlegend = None

        return self.plotGrid(plots,
                             orderlegend=orderlegend,
                             vlines=vlines_axkey,
                             **kwargs,
                             )

    def plotViruses(self,
                    *,
                    ncol=4,
                    nrow=None,
                    sera='all',
                    viruses='all',
                    ignore_virus_serum=None,
                    colors=CBPALETTE,
                    markers=CBMARKERS,
                    serum_to_color_marker=None,
                    max_sera_per_subplot=5,
                    multi_virus_subplots=True,
                    all_subplots=(),
                    titles=None,
                    vlines=None,
                    **kwargs,
                    ):
        """Plot grid with replicate-average of sera for each virus.

        Args:
            `ncol`, `nrow` (int or `None`)
                Specify one of these to set number of columns or rows,
                other should be `None`.
            `sera` ('all' or list)
                Sera to include on plot, in this order, unless one is
                specified in `all_subplots`.
            `viruses` ('all' or list)
                Viruses to include on plot, in this order.
            `ignore_virus_serum` (`None` or dict)
                Specific virus / serum combinations to ignore (not plot). Key
                by virus, and then list sera to ignore.
            `colors` (iterable)
                List of colors for different sera.
            `markers` (iterable)
                List of markers for different sera.
            `serum_to_color_marker` (dict or `None`)
                Optionally specify a specific color and for each serum as
                2-tuples `(color, marker)`. If you use this option, `colors`
                and `markers` are ignored.
            `max_sera_per_subplot` (int)
                Maximum number of sera to show on any subplot.
            `multi_virus_subplots` (bool)
                If a virus has more than `max_sera_per_subplot` sera,
                do we make multiple subplots for it or raise an error?
            `all_subplots` (iterable)
                If making multiple subplots for virus, which sera
                do we show on all subplots? These are also shown first.
            `titles` (`None` or list)
                Specify custom titles for each subplot different than
                `viruses`.
            `vlines` (`None` or dict)
                Add vertical lines to plots. Keyed by virus name, values
                are lists of dicts with a key 'x' giving x-location of vertical
                line, and optional keys 'linewidth', 'color', and 'linestyle'.
            `**kwargs`
                Other keyword arguments that can be passed to
                :meth:`CurveFits.plotGrid`.

        Returns:
            The 2-tuple `(fig, axes)` of matplotlib figure and 2D axes array.

        """
        sera, viruses = self._sera_viruses_lists(sera, viruses)
        viruses = list(collections.OrderedDict.fromkeys(viruses))

        if titles is None:
            titles = viruses
        elif len(viruses) != len(titles):
            raise ValueError(f"`titles`, `viruses` != length:\n"
                             f"{titles}\n{viruses}")

        if max_sera_per_subplot < 1:
            raise ValueError('`max_sera_per_subplot` must be at least 1')

        # get color scheme for sera
        if serum_to_color_marker:
            extra_sera = set(sera) - set(serum_to_color_marker.keys())
            if extra_sera:
                raise ValueError('sera not in `serum_to_color_marker`: ' +
                                 str(extra_sera))
        elif len(sera) <= min(len(colors), len(markers)):
            # can share scheme among subplots
            ordered_sera = ([s for s in sera if s in all_subplots] +
                            [s for s in sera if s not in all_subplots])
            serum_to_color_marker = {s: (c, m) for (s, c, m) in
                                     zip(ordered_sera, colors, markers)}
        elif min(len(colors), len(markers)) < max_sera_per_subplot:
            raise ValueError('`max_sera_per_subplot` larger than '
                             'number of colors or markers')
        else:
            serum_to_color_marker = None

        # Build a list of plots appropriate for `plotGrid`.
        # Code is complicated because we could have several curve
        # per virus, and in that case need to share sera in
        # `all_subplots` among curves.
        virus_sera = {v: [s for s in self.sera if v in self.viruses[s]]
                      for v in self.allviruses}
        plotlist = []
        vlines_list = []
        for virus, title in zip(viruses, titles):
            if ignore_virus_serum and virus in ignore_virus_serum:
                ignore_serum = ignore_virus_serum[virus]
            else:
                ignore_serum = {}
            curvelist = []
            iserum = 0
            virus_shared_sera = [s for s in virus_sera[virus] if
                                 (s in sera) and (s in all_subplots) and
                                 (s not in ignore_serum)]
            virus_unshared_sera = [s for s in virus_sera[virus] if
                                   (s in sera) and
                                   (s not in all_subplots) and
                                   (s not in ignore_serum)]
            unshared = int(bool(len(virus_unshared_sera)))
            if len(virus_shared_sera) > max_sera_per_subplot - unshared:
                raise ValueError(f"virus {virus} has too many subplot-shared "
                                 'sera (in `all_subplots`) relative to '
                                 'value of `max_sera_per_subplot`:\n'
                                 f"{virus_shared_sera} is more than "
                                 f"{max_sera_per_subplot} viruses.")
            shared_curvelist = []
            for serum in virus_shared_sera + virus_unshared_sera:
                if iserum >= max_sera_per_subplot:
                    if multi_virus_subplots:
                        plotlist.append((title, curvelist))
                        if vlines and (virus in vlines):
                            vlines_list.append(vlines[virus])
                        else:
                            vlines_list.append(None)
                        curvelist = list(shared_curvelist)
                        iserum = len(curvelist)
                        assert iserum < max_sera_per_subplot
                    else:
                        raise ValueError(f"virus {virus} has more than "
                                         '`max_sera_per_subplot` viruses '
                                         'and `multi_virus_subplots` is False')
                if serum_to_color_marker:
                    color, marker = serum_to_color_marker[serum]
                else:
                    color = colors[iserum]
                    marker = markers[iserum]
                curvelist.append({'serum': serum,
                                  'virus': virus,
                                  'replicate': 'average',
                                  'label': serum,
                                  'color': color,
                                  'marker': marker,
                                  })
                if serum in virus_shared_sera:
                    shared_curvelist.append(curvelist[-1])
                iserum += 1
            if curvelist:
                plotlist.append((title, curvelist))
                if vlines and (virus in vlines):
                    vlines_list.append(vlines[virus])
                else:
                    vlines_list.append(None)
        if not plotlist:
            raise ValueError('no curves for these viruses / sera')

        # get number of columns
        if (nrow is not None) and (ncol is not None):
            raise ValueError('either `ncol` or `nrow` must be `None`')
        elif isinstance(nrow, int) and nrow > 0:
            ncol = math.ceil(len(plotlist) / nrow)
        elif not (isinstance(ncol, int) and ncol > 0):
            raise ValueError('`nrow` or `ncol` must be integer > 0')

        # convert plotlist to plots dict for `plotGrid`
        plots = {}
        vlines_axkey = {}
        assert len(plotlist) == len(vlines_list)
        for iplot, (plot, ivline) in enumerate(zip(plotlist, vlines_list)):
            irow = iplot // ncol
            icol = iplot % ncol
            plots[(irow, icol)] = plot
            if ivline:
                vlines_axkey[(irow, icol)] = ivline

        if serum_to_color_marker and 'orderlegend' not in kwargs:
            orderlegend = serum_to_color_marker.keys()
        else:
            orderlegend = None

        return self.plotGrid(plots,
                             orderlegend=orderlegend,
                             vlines=vlines_axkey,
                             **kwargs,
                             )

    def plotAverages(self,
                     *,
                     color='black',
                     marker='o',
                     **kwargs,
                     ):
        """Plot grid with a curve for each serum / virus pair.

        Args:
            `color` (str)
                Color the curves.
            `marker` (str)
                Marker for the curves.
            `**kwargs`
                Other keyword arguments that can be passed to
                :meth:`CurveFits.plotReplicates`.

        Returns:
            The 2-tuple `(fig, axes)` of matplotlib figure and 2D axes array.

        """
        return self.plotReplicates(average_only=True,
                                   colors=[color],
                                   markers=[marker],
                                   **kwargs)

    def plotReplicates(self,
                       *,
                       ncol=4,
                       nrow=None,
                       sera='all',
                       viruses='all',
                       colors=CBPALETTE,
                       markers=CBMARKERS,
                       subplot_titles='{serum} vs {virus}',
                       show_average=False,
                       average_only=False,
                       **kwargs,
                       ):
        """Plot grid with replicates for each serum / virus on same plot.

        Args:
            `ncol`, `nrow` (int or `None`)
                Specify one of these to set number of columns or rows.
            `sera` ('all' or list)
                Sera to include on plot, in this order.
            `viruses` ('all' or list)
                Viruses to include on plot, in this order.
            `colors` (iterable)
                List of colors for different replicates.
            `markers` (iterable)
                List of markers for different replicates.
            `subplot_titles` (str)
                Format string to build subplot titles from *serum* and *virus*.
            `show_average` (bool)
                Include the replicate-average as a "replicate" in plots.
            `average_only` (bool)
                Show **only** the replicate-average on each plot. No
                legend in this case.
            `**kwargs`
                Other keyword arguments that can be passed to
                :meth:`CurveFits.plotGrid`.

        Returns:
            The 2-tuple `(fig, axes)` of matplotlib figure and 2D axes array.

        """
        try:
            subplot_titles.format(virus='dummy', serum='dummy')
        except KeyError:
            raise ValueError(f"`subplot_titles` {subplot_titles} invalid. "
                             'Should have format keys only for virus '
                             'and serum')

        sera, viruses = self._sera_viruses_lists(sera, viruses)

        # get replicates and make sure there aren't too many
        nplottable = max(len(colors), len(markers))
        if average_only:
            replicates = ['average']
        else:
            replicates = collections.OrderedDict()
            if show_average:
                replicates['average'] = True
            for serum, virus in itertools.product(sera, viruses):
                if virus in self.viruses[serum]:
                    for replicate in self.replicates[(serum, virus)]:
                        if replicate != 'average':
                            replicates[replicate] = True
            replicates = list(collections.OrderedDict(replicates).keys())
        if len(replicates) > nplottable:
            raise ValueError('Too many unique replicates. There are'
                             f"{len(replicates)} ({', '.join(replicates)}) "
                             f"but only {nplottable} `colors` or `markers`.")

        # build list of plots appropriate for `plotGrid`
        plotlist = []
        for serum, virus in itertools.product(sera, viruses):
            if virus in self.viruses[serum]:
                title = subplot_titles.format(serum=serum, virus=virus)
                curvelist = []
                for i, replicate in enumerate(replicates):
                    if replicate in self.replicates[(serum, virus)]:
                        curvelist.append({'serum': serum,
                                          'virus': virus,
                                          'replicate': replicate,
                                          'label': {False: replicate,
                                                    True: None}[average_only],
                                          'color': colors[i],
                                          'marker': markers[i],
                                          })
                if curvelist:
                    plotlist.append((title, curvelist))
        if not plotlist:
            raise ValueError('no curves for these sera / viruses')

        # get number of columns
        if (nrow is not None) and (ncol is not None):
            raise ValueError('either `ncol` or `nrow` must be `None`')
        elif isinstance(nrow, int) and nrow > 0:
            ncol = math.ceil(len(plotlist) / nrow)
        elif not (isinstance(ncol, int) and ncol > 0):
            raise ValueError('`nrow` or `ncol` must be integer > 0')

        # convert plotlist to plots dict for `plotGrid`
        plots = {}
        for iplot, plot in enumerate(plotlist):
            plots[(iplot // ncol, iplot % ncol)] = plot

        return self.plotGrid(plots, **kwargs)

    def _sera_viruses_lists(self, sera, viruses):
        """Check and build lists of `sera` and their `viruses`.

        Args:
            `sera` ('all' or list)
            `viruses` ('all' or list)

        Returns:
            The 2-tuple `(sera, viruses)` which are checked lists.

        """
        if isinstance(sera, str) and sera == 'all':
            sera = self.sera
        else:
            extra_sera = set(sera) - set(self.sera)
            if extra_sera:
                raise ValueError(f"unrecognized sera: {extra_sera}")

        allviruses = collections.OrderedDict()
        for serum in sera:
            for virus in self.viruses[serum]:
                allviruses[virus] = True
        allviruses = list(allviruses.keys())

        if isinstance(viruses, str) and viruses == 'all':
            viruses = allviruses
        else:
            extra_viruses = set(viruses) - set(allviruses)
            if extra_viruses:
                raise ValueError('unrecognized viruses for specified '
                                 f"sera: {extra_viruses}")

        return sera, viruses

    def plotGrid(self,
                 plots,
                 *,
                 xlabel=None,
                 ylabel=None,
                 widthscale=1,
                 heightscale=1,
                 attempt_shared_legend=True,
                 fix_lims=None,
                 bound_ymin=0,
                 bound_ymax=1,
                 extend_lim=0.07,
                 markersize=6,
                 linewidth=1,
                 linestyle='-',
                 legendtitle=None,
                 orderlegend=None,
                 titlesize=14,
                 labelsize=15,
                 ticksize=12,
                 legendfontsize=12,
                 align_to_dmslogo_facet=False,
                 despine=False,
                 yticklocs=None,
                 sharex=True,
                 sharey=True,
                 vlines=None,
                 ):
        """Plot arbitrary grid of curves.

        Args:
            `plots` (dict)
                Plots to draw on grid. Keyed by 2-tuples `(irow, icol)`, which
                give row and column (0, 1, ... numbering) where plot should be
                drawn. Values are the 2-tuples `(title, curvelist)` where
                `title` is title for this plot (or `None`) and `curvelist`
                is a list of dicts keyed by:

                  - 'serum'
                  - 'virus'
                  - 'replicate'
                  - 'label': label for this curve in legend, or `None`
                  - 'color'
                  - 'marker': https://matplotlib.org/api/markers_api.html

            `xlabel`, `ylabel` (`None`, str, or list)
                Labels for x- and y-axes. If `None`, use `conc_col`
                and `fracinf_col`, respectively. If str, use this shared
                for all axes. If list, should be same length as `plots`
                and gives axis label for each subplot.
            `widthscale`, `heightscale` (float)
                Scale width or height of figure by this much.
            `attempt_shared_legend` (bool)
                Share a single legend among plots if they all share
                in common the same label assigned to the same color / marker.
            `fix_lims` (dict or `None`)
                To fix axis limits, specify any of 'xmin', 'xmax', 'ymin',
                or 'ymax' with specified limit.
            `bound_ymin`, `bound_ymax` (float or `None`)
                Make y-axis min and max at least this small / large.
                Ignored if using `fix_lims` for that axis limit.
            `extend_lim` (float)
                For all axis limits not in `fix_lims`, extend this fraction
                of range above and below bounds / data limits.
            `markersize` (float)
                Size of point marker.
            `linewidth` (float)
                Width of line.
            `linestyle` (str)
                Line style.
            `legendtitle` (str or `None`)
                Title of legend.
            `orderlegend` (`None` or list)
                If specified, place legend labels in this order.
            `titlesize` (float)
                Size of subplot title font.
            `labelsize` (float)
                Size of axis label font.
            `ticksize` (float)
                Size of axis tick fonts.
            `legendfontsize` (float)
                Size of legend fonts.
            `align_to_dmslogo_facet` (`False` or dict)
                Make plot vertically alignable to ``dmslogo.facet_plot``
                with same number of rows; dict should have keys for
                `height_per_ax`, `hspace`, `tmargin`, and `bmargin` with
                same meaning as ``dmslogo.facet_plot``. Also
                `right` and `left` for passing to ``subplots_adjust``.
            `despine` (bool)
                Remove top and right spines from plots.
            `yticklocs` (`None` or list)
                Same meaning as for :meth:`neutcurve.hillcurve.HillCurve.plot`.
            `sharex` (bool)
                Share x-axis scale among plots.
            `sharey` (bool)
                Share y-axis scale among plots.
            `vlines` (dict or `None`)
                Vertical lines to draw. Keyed by 2-tuples `(irow, icol)`, which
                give row and column of plot in grid (0, 1, ... numbering).
                Values are lists of dicts with a key 'x' giving the x-location
                of the vertical line, and optionally keys 'linewidth',
                'color', and 'linestyle'.

        Returns:
            The 2-tuple `(fig, axes)` of matplotlib figure and 2D axes array.

        """
        vline_defaults = {'linewidth': 1.5,
                          'color': 'gray',
                          'linestyle': '--',
                          }

        if not plots:
            raise ValueError('empty `plots`')

        # get number of rows / cols, curves, and data limits
        nrows = ncols = None
        if fix_lims is None:
            fix_lims = {}
        lims = {key: {} for key in plots.keys()}
        for (irow, icol), (_title, curvelist) in plots.items():
            if irow < 0:
                raise ValueError('invalid row index `irow` < 0')
            if icol < 0:
                raise ValueError('invalid row index `icol` < 0')
            if nrows is None:
                nrows = irow + 1
            else:
                nrows = max(nrows, irow + 1)
            if ncols is None:
                ncols = icol + 1
            else:
                ncols = max(ncols, icol + 1)
            for curvedict in curvelist:
                curve = self.getCurve(serum=curvedict['serum'],
                                      virus=curvedict['virus'],
                                      replicate=curvedict['replicate']
                                      )
                curvedict['curve'] = curve
                for lim, attr, f in [('xmin', 'cs', min), ('xmax', 'cs', max),
                                     ('ymin', 'fs', min), ('ymax', 'fs', max)]:
                    if lim in fix_lims:
                        lims[(irow, icol)][lim] = fix_lims[lim]
                    else:
                        val = f(getattr(curve, attr))
                        if lim in lims[(irow, icol)]:
                            val = f(val, lims[(irow, icol)][lim])
                        if lim == 'ymin' and (bound_ymin is not None):
                            lims[(irow, icol)][lim] = min(val, bound_ymin)
                        elif lim == 'ymax' and (bound_ymax is not None):
                            lims[(irow, icol)][lim] = max(val, bound_ymax)
                        else:
                            lims[(irow, icol)][lim] = val

        for share, axtype in [(sharex, 'x'), (sharey, 'y')]:
            if share:
                for limtype, limfunc in [('min', min), ('max', max)]:
                    lim = limfunc(lims[key][axtype + limtype] for key in lims)
                    for key in lims.keys():
                        lims[key][axtype + limtype] = lim

        # check and then extend limits
        for key in plots.keys():
            if lims[key]['xmin'] <= 0:
                raise ValueError('xmin <= 0, which is not allowed')
            yextent = lims[key]['ymax'] - lims[key]['ymin']
            if yextent <= 0:
                raise ValueError('no positive extent for y-axis')
            if 'ymin' not in fix_lims:
                lims[key]['ymin'] -= yextent * extend_lim
            if 'ymax' not in fix_lims:
                lims[key]['ymax'] += yextent * extend_lim
            xextent = math.log(lims[key]['xmax']) - math.log(lims[key]['xmin'])
            if xextent <= 0:
                raise ValueError('no positive extent for x-axis')
            if 'xmin' not in fix_lims:
                lims[key]['xmin'] = math.exp(math.log(lims[key]['xmin']) -
                                             xextent * extend_lim)
            if 'xmax' not in fix_lims:
                lims[key]['xmax'] = math.exp(math.log(lims[key]['xmax']) +
                                             xextent * extend_lim)

        if align_to_dmslogo_facet:
            hparams = dmslogo.facet.height_params(
                                nrows,
                                align_to_dmslogo_facet['height_per_ax'],
                                align_to_dmslogo_facet['hspace'],
                                align_to_dmslogo_facet['tmargin'],
                                align_to_dmslogo_facet['bmargin'],
                                )
            height = hparams['height']
        else:
            height = (1 + 2.25 * nrows) * heightscale

        width = (1 + 3 * ncols) * widthscale
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 sharex=sharex,
                                 sharey=sharey,
                                 squeeze=False,
                                 figsize=(width, height),
                                 )

        # set limits for share axes
        for irow, icol in plots.keys():
            axes[irow, icol].set_xlim(lims[irow, icol]['xmin'],
                                      lims[irow, icol]['xmax'])
            axes[irow, icol].set_ylim(lims[irow, icol]['ymin'],
                                      lims[irow, icol]['ymax'])

        # make plots
        shared_legend = attempt_shared_legend
        kwargs_tup_to_label = {}  # used to determine if shared legend
        legend_handles = collections.defaultdict(list)
        shared_legend_handles = []  # handles if using shared legend
        for i, ((irow, icol), (title, curvelist)) in enumerate(plots.items()):
            ax = axes[irow, icol]
            ax.set_title(title, fontsize=titlesize)
            for curvedict in curvelist:
                kwargs = {'color': curvedict['color'],
                          'marker': curvedict['marker'],
                          'linestyle': linestyle,
                          'linewidth': linewidth,
                          'markersize': markersize,
                          }
                if isinstance(xlabel, list):
                    ixlabel = xlabel[i]
                else:
                    ixlabel = None
                if isinstance(ylabel, list):
                    iylabel = ylabel[i]
                else:
                    iylabel = None
                curvedict['curve'].plot(ax=ax,
                                        xlabel=ixlabel,
                                        ylabel=iylabel,
                                        yticklocs=yticklocs,
                                        **kwargs,
                                        )
                label = curvedict['label']
                if label:
                    handle = Line2D(xdata=[],
                                    ydata=[],
                                    label=label,
                                    **kwargs,
                                    )
                    legend_handles[(irow, icol)].append(handle)
                    if shared_legend:
                        kwargs_tup = tuple(sorted(kwargs.items()))
                        if kwargs_tup in kwargs_tup_to_label:
                            if kwargs_tup_to_label[kwargs_tup] != label:
                                shared_legend = False
                        else:
                            kwargs_tup_to_label[kwargs_tup] = label
                            shared_legend_handles.append(handle)
            ax.tick_params('both', labelsize=ticksize, bottom=True, left=True,
                           right=False, top=False)
            if despine:
                dmslogo.utils.despine(ax=ax)
            if vlines and ((irow, icol) in vlines):
                for vline in vlines[(irow, icol)]:
                    vline_d = vline_defaults.copy()
                    for key, val in vline.items():
                        vline_d[key] = val
                    ax.axvline(vline_d['x'],
                               linestyle=vline_d['linestyle'],
                               linewidth=vline_d['linewidth'],
                               color=vline_d['color'])

        # draw legend(s)
        legend_kwargs = {'fontsize': legendfontsize,
                         'numpoints': 1,
                         'markerscale': 1,
                         'handlelength': 1,
                         'labelspacing': 0.1,
                         'handletextpad': 0.4,
                         'frameon': True,
                         'borderaxespad': 0.1,
                         'borderpad': 0.2,
                         'title': legendtitle,
                         'title_fontsize': legendfontsize + 1,
                         'framealpha': 0.6,
                         }

        def _ordered_legend(hs):
            """Get ordered legend handles."""
            if not orderlegend:
                return hs
            else:
                order_dict = {h: i for i, h in enumerate(orderlegend)}
                h_labels = [h.get_label() for h in hs]
                extra_hs = set(h_labels) - set(orderlegend)
                if extra_hs:
                    raise ValueError('there are legend handles not in '
                                     f"`orderlegend`: {extra_hs}")
                return [h for _, h in sorted(zip(h_labels, hs),
                                             key=lambda x: order_dict[x[0]])]

        if shared_legend and shared_legend_handles:
            if align_to_dmslogo_facet:
                right = align_to_dmslogo_facet['right']
                ranchor = right + 0.15 * (1 - right)
            else:
                ranchor = 1
            shared_legend_handles = _ordered_legend(shared_legend_handles)
            # shared legend as here: https://stackoverflow.com/a/17328230
            fig.legend(handles=shared_legend_handles,
                       labels=[h.get_label() for h in shared_legend_handles],
                       loc='center left',
                       bbox_to_anchor=(ranchor, 0.5),
                       bbox_transform=fig.transFigure,
                       **legend_kwargs,
                       )
        elif legend_handles:
            for (irow, icol), handles in legend_handles.items():
                ax = axes[irow, icol]
                handles = _ordered_legend(handles)
                ax.legend(handles=handles,
                          labels=[h.get_label() for h in handles],
                          loc='lower left',
                          **legend_kwargs,
                          )

        # hide unused axes
        for irow, icol in itertools.product(range(nrows), range(ncols)):
            if (irow, icol) not in plots:
                axes[irow, icol].set_axis_off()

        # common axis labels as here: https://stackoverflow.com/a/53172335
        bigax = fig.add_subplot(111, frameon=False)
        bigax.grid(False)
        bigax.tick_params(labelcolor='none', top=False, bottom=False,
                          left=False, right=False, which='both')
        if xlabel is None:
            bigax.set_xlabel(self.conc_col, fontsize=labelsize, labelpad=10)
        elif not isinstance(xlabel, list):
            bigax.set_xlabel(xlabel, fontsize=labelsize, labelpad=10)
        if ylabel is None:
            bigax.set_ylabel(self.fracinf_col, fontsize=labelsize, labelpad=10)
        elif not isinstance(ylabel, list):
            bigax.set_ylabel(ylabel, fontsize=labelsize, labelpad=10)

        if align_to_dmslogo_facet:
            fig.subplots_adjust(hspace=hparams['hspace'],
                                top=hparams['top'],
                                bottom=hparams['bottom'],
                                left=align_to_dmslogo_facet['left'],
                                right=align_to_dmslogo_facet['right'],
                                )
        else:
            fig.tight_layout()

        return fig, axes


if __name__ == '__main__':
    import doctest
    doctest.testmod()
