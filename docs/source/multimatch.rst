multimatch
==========


The ``multimatch`` command is the standalone equivalent of the MultiMatch
toolbox and is easiest executed directly from the command line.

execution via commandline
^^^^^^^^^^^^^^^^^^^^^^^^^

The computation of the similarity between two scanpaths doesn't involve anything
beyond the command line keyword ``multimatch`` followed by two inputs,
corresponding to tabseparated files with a fixation vector:


.. code::

   multimatch path/to/scanpath_one path/to/scanpath_two

``multimatch`` needs the screensize of the presentation screen as an input. The
default is 1280 x 720px. If this needs adjustment, use the optional --screensize
flag:
- ``--screensize``: in px. Specify first x, then y dimension. Default: 1280 x
720.

.. code::

   multimatch path/to/scanpath_one path/to/scanpath_two --screensize 1280 720


Optionally, scanpaths can be simplified to reduce their complexity. To simplify
scanpaths, specify the following arguments:

- ``--direction_threshold``: If two consecutive saccades have a small angle, they will be
  combined. Should be in degrees, such as ``45.0`` for 45°
- ``--amplitude_threshold``: If two consecutive saccades are short, they will be
  combines. Should be in pixel, such as ``100.0`` for 100px.
- ``--duration_threshold``: Only if the intermediate fixation's durations are
  shorter than this threshold the above simplification will be performed. Should
  be in seconds, such as ``0.1`` for 100ms.

**Note**: If either direction- or amplitude threshold is specified as 0, no
grouping will be performed!


A commandline call of the module **with** simplification would hence look like
this:

.. code::

   multimatch path/to/scanpath_one path/to/scanpath_two --screensize 1280 720
   --direction_threshold 45.0 --amplitude_threshold 100.0 --direction_threshold 0.1


There are no guidelines how much simplification is appropriate, and it is strongly dependent
on individual use case. The original matlab toolbox implements a default
amplitude threshold of 10% of the screen diagonal as amplitude, 45° as angle, and 300ms as
duation thresholds.


execution within a python instance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you wish to use the functionality of multimatch within a running python
instance such as ipython, you can import the module and use the function
``docomparison``. Here is an example:

.. code::

   import multimatch as m
   import numpy as np

   #read in data
   fix_vector1 = np.recfromcsv('data/fixvectors/segment_0_sub-01.tsv',
   delimiter='\t', dtype={'names': ('start_x', 'start_y', 'duration'),
   'formats': ('f8', 'f8', 'f8')})
   fix_vector2 = np.recfromcsv('data/fixvectors/segment_0_sub-19.tsv',
   delimiter='\t', dtype={'names': ('start_x', 'start_y', 'duration'),
   'formats': ('f8', 'f8', 'f8')})

   #execution with multimatch's docomparison() function without grouping
   m.docomparison(fix_vector1, fix_vector2, sz=[1280, 720])

   #execution with multimatch's docomparison() function with grouping
   m.docomparison(fix_vector1, fix_vector2, sz=[1280, 720], grouping=True, TDir=30.0,
   TDur=0.1, TAmp=100.1)



..  TODO .. automodule:: multimatch.multimatch
    :members:
    :undoc-members:
    :show-inheritance:
..  TODO \n .. automodule:: multimatch.multimatch

