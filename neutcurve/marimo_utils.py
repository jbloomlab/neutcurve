"""
============
marimo_utils
============

Utilities for working with `marimo <https://marimo.io/>`_ notebooks.
"""

import base64
import io

import matplotlib
import matplotlib.figure


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

    import marimo as mo

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
        return mo.Html(
            f"""
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
        )

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
        buf.seek(0)
        return mo.pdf(src=buf, width="100%", height="80vh")

    elif display_method == "png8":
        import PIL

        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=80, metadata={})
        im = PIL.Image.open(io.BytesIO(buf.getvalue())).quantize(
            colors=48, dither=PIL.Image.NONE
        )
        out = io.BytesIO()
        im.save(out, format="PNG", optimize=True, pnginfo=PIL.PngImagePlugin.PngInfo())
        data = base64.b64encode(out.getvalue()).decode("ascii")
        return mo.Html(f'<img src="data:image/png;base64,{data}" alt="figure">')

    else:
        raise ValueError(f"Invalid {display_method=}")


if __name__ == "__main__":
    import doctest

    doctest.testmod()
