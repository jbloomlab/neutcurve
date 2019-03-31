"""
==================
curvefits
==================
Defines :class:`CurveFits` to fit curves and display / plot results.
"""


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
        `replicate_col` (str or `None`)
            If there are replicate measurements for serum / virus pairs,
            should give column indicating the replicate. Replicates must all
            have the same concentrations for each serum / virus combination.
        `fixbottom` (`False` or float)
            Same meaning as for :mod:`hillcurve.HillCurve`.
        `fixtop` (`False` or float)
            Same meaning as for :mod:`hillcurve.HillCurve`.

    Attributes of a :class:`CurveFits` include all args plus:
        `sera` (list)
            List of all serum names in `serum_col` of `data`, in order
            they occur in `data`.
        `viruses` (dict)
            For each serum in `sera`, `viruses[serum]` gives all viruses
            for that serum in the order they occur in `data`.
        `replicates` (`None` or dict)
            If `replicate_col` is not `None`, then
            `replicates[(serum, virus)]` is list of all replicates for
            that serum and virus in the order they occur in data.
    """

    def __init__(self,
                 data,
                 *,
                 conc_col='concentration',
                 fracinf_col='fraction infectivity',
                 serum_col='serum',
                 virus_col='virus',
                 replicate_col=None,
                 fixbottom=False,
                 fixtop=1,
                 ):
        """See main class docstring."""
        # make args into attributes
        self.data = data
        self.conc_col = conc_col
        self.fracinf_col = fracinf_col
        self.serum_col = serum_col
        self.virus_col = virus_col
        self.replicate_col = replicate_col
        self.fixbottom = fixbottom
        self.fixtop = fixtop

        # check for required columns
        cols = [self.conc_col, self.fracinf_col, self.serum_col,
                self.virus_col]
        if self.replicate_col is not None:
            cols.append(self.replicate_col)
        if len(cols) != len(set(cols)):
            raise ValueError('duplicate column names:\n\t' + '\n\t'.join(cols))
        if not (set(cols) <= set(data.columns)):
            raise ValueError('`data` lacks required columns, which are:\n\t' +
                             '\n\t'.join(cols))

        # create sera / viruses / replicates attributes
        self.sera = self.data[self.serum_col].unique().tolist()
        self.viruses = {serum: (data
                                .query(f"{self.serum_col} == {serum}")
                                [self.virus_col]
                                .unique()
                                .tolist()
                                )
                        for serum in self.sera}
        if self.replicate_col is None:
            self.replicates = None
        else:
            self.replicates = {(s, v): (data
                                        .query(f"({self.serum_col} == {s}) & "
                                               f"({self.virus_col} == {v})")
                                        [self.replicate_col]
                                        .unique()
                                        .tolist()
                                        )
                               for s, v in self.viruses.items()}


if __name__ == '__main__':
    import doctest
    doctest.testmod()
