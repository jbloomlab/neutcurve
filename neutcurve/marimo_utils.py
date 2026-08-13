"""
============
marimo_utils
============

Utilities for working with `marimo <https://marimo.io/>`_ notebooks.
"""

import matplotlib
import matplotlib.figure

from neutcurve.fig_utils import fig_html


def display_fig_marimo(fig, display_method):
    """Display large matplotlib figure via marimo.

    This function is designed when you have a very large matplotlib figure (eg, as
    produced by the plotting functions of :class:`neutcurve.curvefits.CurveFits`)
    and you want to display it in marimo notebook.

    Running this function requires you to have separately installed
    `marimo <https://marimo.io/>`_ and (if you are using `display_method="png8")
    `pillow <https://python-pillow.github.io/>`_

    Args:
        `fig` (matplotlib..figure.Figure)
            The figure we want to display.
        `display_method` {"inline", "svg", "pdf", "png8"}
            Display the figure just inline, as a SVG, as a PDF, or as a PNG8.
            In general, displaying as a PNG8 will be the smallest size although
            also the lowest resolution.

    Returns:
        output_obj
            The returned object can be display in marimo, via
            ``marimo.output.append(output_obj)``

    """
    if not isinstance(fig, matplotlib.figure.Figure):
        raise ValueError(
            f"Expected `fig` to be matplotlib.figure.Figure, instead {type(fig)=}"
        )

    if display_method == "inline":
        return fig

    if display_method not in {"svg", "pdf", "png8"}:
        raise ValueError(
            f"Invalid {display_method=}, valid methods are "
            "'inline', 'svg', 'pdf', 'png8'"
        )

    import marimo as mo

    return mo.Html(fig_html(fig, display_method))


if __name__ == "__main__":
    import doctest

    doctest.testmod()
