"""
============
colorschemes
============

Color schemes.
"""

#: color-blind safe palette with gray, from
#: http://bconnelly.net/2013/10/creating-colorblind-friendly-figures
CBPALETTE = ('#999999', '#E69F00', '#56B4E9', '#009E73',
             '#F0E442', '#0072B2', '#D55E00', '#CC79A7')

#: color-blind safe palette with black, from
#: http://bconnelly.net/2013/10/creating-colorblind-friendly-figures
CBBPALETTE = ('#000000', '#E69F00', '#56B4E9', '#009E73',
              '#F0E442', '#0072B2', '#D55E00', '#CC79A7')

#: list of matplotlib markers same length as `CBPALETTE`, from
#: https://matplotlib.org/api/markers_api.html
CBMARKERS = ('o', '^', 's', 'D', 'v', '<', '>', 'p')
assert len(CBPALETTE) == len(CBBPALETTE) == len(CBMARKERS)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
