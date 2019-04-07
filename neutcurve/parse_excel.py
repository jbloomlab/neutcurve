"""
==================
parse_excel
==================

Defines functions for parsing data directly from Excel in specific formats.
"""

import os

import pandas as pd


def parseRachelStyle2019(*,
                         excelfile,
                         sheet_mapping,
                         ):
    """Parse data from Rachel Eguia's 2019-style assays on Bloom lab Tecan.

    Args:
        `excelfile` (str)
            Excel file with data exactly as produced by plate reader.
        `sheet_mapping` (dict)
            Describes data in each sheet of `excelfile`. There should
            be a key matching the name of each Excel sheet. The value
            for each key is another dict with the following keys / values
            describing the data for that plate (we assume one sheet per
            serum  / virus pair):

              - 'serum': name of serum.
              - 'virus': name of virus.
              - 'dilution_factor': fold-dilution of serum along plate.
              - 'initial_concentration': final serum concentration in
                the first well used to start the dilution series.

    Returns:
        A pandas DataFrame in the tidy format read by
        :class:`neutcurve.curvefits.CurveFits`.

    """
    if not os.path.isfile(excelfile):
        raise ValueError(f"cannot find `excelfile` {excelfile}")

    sheet_data = pd.read_excel(excelfile,
                               sheet_name=None,  # read all sheets
                               engine='xlrd',
                               skiprows=range(0, 30),
                               index_col=0,
                               nrows=8,
                               )

    # keys of sheet_mapping should be str
    sheet_mapping = {str(key): val for key, val in sheet_mapping.items()}

    extra_sheets = set(sheet_data) - set(sheet_mapping)
    if extra_sheets:
        raise ValueError(f"`excelfile` {excelfile} has the following extra "
                         f"sheets not in `sheet_mapping`: {extra_sheets}")

    neutdata = []
    for sheet, sheet_info in sheet_mapping.items():

        if sheet not in sheet_data:
            raise ValueError(f"sheet {sheet} specified in `sheet_mapping` "
                             f"is not present in `excelfile` {excelfile}")

        kwargs = {'layout': 'RachelStyle2019'}
        for key in ['dilution_factor', 'initial_concentration',
                    'serum', 'virus']:
            try:
                kwargs[key] = sheet_info[key]
            except KeyError:
                raise ValueError(f"`sheet_mapping` does not specify {key} "
                                 f"for sheet {sheet}")

        neutdata.append(_parse_BloomLab_TecanPlate(
                sheet_df=sheet_data[sheet],
                sheet_name_for_err=f"sheet {sheet} of {excelfile}",
                **kwargs,
                ))

    return pd.concat(neutdata, ignore_index=True, sort=False)


def _parse_BloomLab_TecanPlate(*,
                               sheet_df,
                               sheet_name_for_err,
                               initial_concentration,
                               dilution_factor,
                               layout,
                               serum,
                               virus,
                               ):
    """Process data frame for one 96-well plate of Bloom lab Tecan."""
    # expected rows / columns
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    cols = list(range(1, 13))

    # check for expected rows / columns
    if sheet_df.index.tolist() != rows:
        raise ValueError(f"{sheet_name_for_err} does not have expected rows.\n"
                         f"Expected: {rows}\nGot: {sheet_df.index.tolist()}")
    if sheet_df.columns.tolist() != cols:
        raise ValueError(f"{sheet_name_for_err} does not have expected "
                         f"columns.\nExpected: {cols}\n"
                         f"Got: {sheet_df.columns.tolist()}")
    if not all(pd.api.types.is_numeric_dtype(sheet_df[col]) for col in cols):
        raise ValueError(f"Entries in {sheet_name_for_err} not all numerical")

    if layout == 'RachelStyle2019':
        reps = {rep: sheet_df.loc[row] for rep, row in [('1', 'D'),
                                                        ('2', 'E'),
                                                        ('3', 'F')]}
        no_serum = (sheet_df.loc['C'] + sheet_df.loc['G']) / 2
        virus_only = sheet_df.loc['B']
        media_only = (sheet_df.loc['A'] + sheet_df.loc['H']) / 2
        if any(media_only > virus_only):
            raise ValueError(f"In {sheet_name_for_err}, the media-only "
                             f"control has more signal than the virus-only "
                             f"control, which is unexpected:\n{sheet_df}")
        # compute fraction infectivity for each replicate
        data = {rep: (repval - virus_only) / (no_serum - virus_only)
                for rep, repval in reps.items()}
        data['concentration'] = [initial_concentration / dilution_factor**i
                                 for i in range(len(cols) - 1, -1, -1)]
        df = (pd.DataFrame(data)
              .melt(id_vars='concentration',
                    var_name='replicate',
                    value_name='fraction infectivity'
                    )
              )
    else:
        raise ValueError(f"invalid `layout` of {layout}")

    return (df
            .assign(serum=serum,
                    virus=virus)
            [['serum', 'virus', 'replicate', 'concentration',
              'fraction infectivity']]
            )


if __name__ == '__main__':
    import doctest
    doctest.testmod()
