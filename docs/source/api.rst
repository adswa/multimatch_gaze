API
===


The ``multimatch`` command is the standalone equivalent of the MultiMatch
toolbox and is easiest executed directly from the command line.


Command line
^^^^^^^^^^^^

The computation of the similarity between two scanpaths doesn't involve anything
beyond the command line keyword ``multimatch`` followed by two input files,
corresponding to tab-separated files with a fixation vector, and the screensize in
pixel, supplied as two consecutive integers corresponding to the x and y dimension
of the screen:


.. code::

   multimatch path/to/scanpath_one path/to/scanpath_two x_dim y_dim

The input files will be read with ``numpy``\s ``recfromcsv()`` and should contain
one fixation per line. The first column is assumed to contain (center) x-coordinates
of the fixation in px, the second column is assumed to contain (center) y-coordinates
in px, and the third column is assumed to contain the fixation duration in seconds.
Example files with this structure can be found here_.

 .. _here: https://github.com/adswa/multimatch/tree/master/data/fixvectors

An examplary command line call that you could execute if you have cloned the
repository looks like this:

.. code::

   multimatch data/fixvectors/segment_0_sub-01.tsv data/fixvectors/segment_0_sub-19.tsv 1280 720

**Scanpath simplification**

Optionally, scanpaths can be simplified to reduce their complexity. To simplify
scanpaths, specify the following arguments:

- ``--direction-threshold``: If two consecutive saccades have a small angle, they will be
  combined. Should be in degrees, such as ``45.0`` for 45°
- ``--amplitude-threshold``: If two consecutive saccades are short, they will be
  combined. Should be in pixel, such as ``100.0`` for 100px.
- ``--duration-threshold``: Only if the intermediate fixation's durations are
  shorter than this threshold the above simplification will be performed. Should
  be in seconds, such as ``0.1`` for 100ms.

**Note**: If either direction- or amplitude threshold is specified as 0, no
grouping will be performed!


A commandline call of the module **with** simplification would hence look like
this:

.. code::

   multimatch data/fixvectors/segment_0_sub-01.tsv data/fixvectors/segment_0_sub-19.tsv 1280 720
   --direction-threshold 45.0 --amplitude-threshold 100.0 --duration-threshold 0.1


There are no guidelines whether and if so, how much,
simplification is appropriate, and it is strongly dependent
on individual use case. The original matlab toolbox implements a default
amplitude threshold of 10% of the screen diagonal as amplitude, 45° as angle, and 300ms as
duration thresholds. ``multimatch`` has defaults of 0 for simplification parameters
(i.e. simplification is not performed by default).

**Output configuration**

The way results are displayed in the command line can be configured with the ``-o``/``--output-type``
parameter.
Three different formats are possible:

- ``hr`` (default): Results are returned row-wise, with dimension name. This is the
  most human readable format, and good for a quick glance at results:

.. code::

   Vector similarity = <value>
   Direction similarity = <value>
   ...

- ``single-row``: Results are returned in a single row, delimited with tabs, and without
  dimension name. Makes it easy to collate results in a table:

.. code::

   <vectorsim>\t<directionsim>\t<lengthsim>\t<positionsim>\t<durationsim>

- ``single-del``: Results are returned row-wise, with tabs seperating dimension name
  and value. This makes it easy to pick out a selection of scores:

.. code::

   vector\t<value>
   direction\t<value>
   length\t<value>
   position\t<value>
   duration\t<value>

Python
^^^^^^

If you wish to use the functionality of multimatch within a running python
instance such as ipython, you can import the module and use the function
``docomparison``. Here is an example:

.. code::

   import multimatch as m
   import numpy as np

   # read in data
   fix_vector1 = np.recfromcsv('data/fixvectors/segment_0_sub-01.tsv',
   delimiter='\t', dtype={'names': ('start_x', 'start_y', 'duration'),
   'formats': ('f8', 'f8', 'f8')})
   fix_vector2 = np.recfromcsv('data/fixvectors/segment_0_sub-19.tsv',
   delimiter='\t', dtype={'names': ('start_x', 'start_y', 'duration'),
   'formats': ('f8', 'f8', 'f8')})

   # execution with multimatch's docomparison() function without grouping
   m.docomparison(fix_vector1, fix_vector2, screensize=[1280, 720])

   # execution with multimatch's docomparison() function with grouping
   m.docomparison(fix_vector1, fix_vector2, screensize=[1280, 720], grouping=True, TDir=30.0,
   TDur=0.1, TAmp=100.1)

The results will be returned as an array, such as ``[0.98, 0.87, 0.45, 0.78, 0.80]``.