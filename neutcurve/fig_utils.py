"""
=========
fig_utils
=========

Utilities for rendering `matplotlib <https://matplotlib.org/>`_ figures as HTML.

"""

import base64
import io

import matplotlib
import matplotlib.figure


def fig_html(fig, display_method):
    """Render a large matplotlib figure as a self-contained HTML string.

    This function is designed for when you have a very large matplotlib figure (eg, as
    produced by the plotting functions of :class:`neutcurve.curvefits.CurveFits`) and you
    want to put it in an HTML page. For ``display_method="svg"`` and ``"pdf"`` the figure
    is placed in a fixed-height scrollable box, so that a figure with many panels does
    not force the rest of the page arbitrarily far down.

    Requires `pillow <https://python-pillow.github.io/>`_ if you are using
    ``display_method="png8"``; it is installed with `matplotlib` in any case.

    Args:
        `fig` (matplotlib.figure.Figure)
            The figure we want to render.
        `display_method` {"svg", "pdf", "png8"}
            Render the figure as a SVG, as a PDF, or as a PNG8. In general, rendering
            as a PNG8 will be the smallest size although also the lowest resolution.

    Returns:
        `str`
            HTML that can be embedded in a page, or passed to ``marimo.Html``.

    Example:

    >>> import matplotlib.figure
    >>> fig = matplotlib.figure.Figure(figsize=(2, 2))
    >>> _ = fig.subplots().plot([0, 1], [0, 1])

    The SVG is wrapped in the scrollable box, and is a real SVG rather than an image:

    >>> html = fig_html(fig, "svg")
    >>> '<div id="svgwrap" style="width:100%;height:80vh;overflow:auto">' in html
    True
    >>> "<svg" in html
    True

    The other methods embed the figure as a data URI:

    >>> fig_html(fig, "pdf").startswith("<iframe src='data:application/pdf;base64,")
    True
    >>> fig_html(fig, "png8").startswith('<img src="data:image/png;base64,')
    True

    Anything else is an error. Note in particular that there is no ``"inline"`` method
    here, as a bare figure has no HTML form; see
    :func:`neutcurve.marimo_utils.display_fig_marimo` for that:

    >>> fig_html(fig, "inline")
    Traceback (most recent call last):
        ...
    ValueError: Invalid display_method='inline', valid methods are 'svg', 'pdf', 'png8'

    """
    if not isinstance(fig, matplotlib.figure.Figure):
        raise ValueError(
            f"Expected `fig` to be matplotlib.figure.Figure, instead {type(fig)=}"
        )

    if display_method == "svg":
        buf = io.BytesIO()
        with matplotlib.rc_context(
            {
                "svg.fonttype": "none",  # keep text as text, not paths
                "svg.image_inline": True,  # embed small images if present
                "svg.hashsalt": "fixed-1",  # deterministic ids in the SVG
                "path.simplify": True,
                "path.simplify_threshold": 0.2,
            }
        ):
            fig.savefig(buf, format="svg", metadata={})
        svg_text = buf.getvalue().decode("utf-8")
        return f"""
<style>
#svgwrap svg {{
  width: 100% !important;
  height: auto !important;
  max-width: 100%;
  display: block;
}}
</style>
<div id="svgwrap" style="width:100%;height:80vh;overflow:auto">
  {svg_text}
</div>
"""

    elif display_method == "pdf":
        buf = io.BytesIO()
        with matplotlib.rc_context(
            {
                "pdf.fonttype": 42,
                "pdf.compression": 7,
                "path.simplify": True,
                "path.simplify_threshold": 0.2,
            }
        ):
            fig.savefig(buf, format="pdf", metadata={})
        data = base64.b64encode(buf.getvalue()).decode("ascii")
        return (
            f"<iframe src='data:application/pdf;base64,{data}' "
            "style='border-radius: 4px;width: 100%;height: 80vh' "
            "frameborder='0'></iframe>"
        )

    elif display_method == "png8":
        import PIL.Image
        import PIL.PngImagePlugin

        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=80, metadata={})
        im = PIL.Image.open(io.BytesIO(buf.getvalue())).quantize(
            colors=48, dither=PIL.Image.Dither.NONE
        )
        out = io.BytesIO()
        im.save(out, format="PNG", optimize=True, pnginfo=PIL.PngImagePlugin.PngInfo())
        data = base64.b64encode(out.getvalue()).decode("ascii")
        return f'<img src="data:image/png;base64,{data}" alt="figure">'

    else:
        raise ValueError(
            f"Invalid {display_method=}, valid methods are 'svg', 'pdf', 'png8'"
        )


if __name__ == "__main__":
    import doctest

    doctest.testmod()
