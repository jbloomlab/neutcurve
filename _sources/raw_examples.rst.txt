.. _raw_examples:

Raw data from plate reader
==================================

``neutcurve`` also contains functions that allow the raw data (i.e., Excel files)
to be directly parsed for some common experimental formats. This parsing
automates creation of the tidy data frame that is passed
to :class:`neutcurve.curvefits.CurveFits` as in :ref:`curvefits_example`.

Here are formats that can be parsed in this way:

.. toctree::
   :maxdepth: 1

   rachelstyle2019_example
