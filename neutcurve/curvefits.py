"""
==================
curvefits
==================
Defines :class:`CurveFits` to fit curves and display / plot results.
"""


import neutcurve


class CurveFits:
    """Fit and display :mod:`hillcurve.HillCurve` curves.

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
            Same meaning as for :mod:`hillcurve.HillCurve`.
        `fixtop` (`False` or float)
            Same meaning as for :mod:`hillcurve.HillCurve`.

    Attributes of a :class:`CurveFits` include all args except `data` plus:
        `df` (pandas DataFrame)
            Copy of `data` that only has relevant columns, has additional rows
            with `replicate_col` of 'average' that hold replicate averages, and
            an additional column 'stdev' holding the standard deviation of
            fraction infectivity for 'average' if there multiple replicates,
            and `nan` for all other replicates.
        `sera` (list)
            List of all serum names in `serum_col` of `data`, in order
            they occur in `data`.
        `viruses` (dict)
            For each serum in `sera`, `viruses[serum]` gives all viruses
            for that serum in the order they occur in `data`.
        `replicates` (dict)
            `replicates[(serum, virus)]` is list of all replicates for
            that serum and virus in the order they occur in `data`.
    """

    def __init__(self,
                 data,
                 *,
                 conc_col='concentration',
                 fracinf_col='fraction infectivity',
                 serum_col='serum',
                 virus_col='virus',
                 replicate_col='replicate',
                 fixbottom=False,
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

        # check for required columns
        cols = [self.conc_col, self.fracinf_col, self.serum_col,
                self.virus_col, self.replicate_col]
        if len(cols) != len(set(cols)):
            raise ValueError('duplicate column names:\n\t' + '\n\t'.join(cols))
        if not (set(cols) <= set(data.columns)):
            raise ValueError('`data` lacks required columns, which are:\n\t' +
                             '\n\t'.join(cols))

        self.df = data[cols]

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
                self.replicates[(serum, virus)] = virus_reps
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

        self._hillcurves = {}  # curves computed by `getHillCurve` go here

    def getHillCurve(self, *, serum, virus, replicate):
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
            A :mod:class:`HillCurve`.

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

            curve = neutcurve.HillCurve(cs=idata[self.conc_col],
                                        fs=idata[self.fracinf_col],
                                        fixbottom=self.fixbottom,
                                        fixtop=self.fixtop,
                                        )

            self._hillcurves[key] = curve

        return self._hillcurves[key]


if __name__ == '__main__':
    import doctest
    doctest.testmod()
