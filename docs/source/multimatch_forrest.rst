multimatch_forrest
==================

The ``multimatch_forrest`` command is additional functionality intended to use
inputs from the studyforrest_ phase 2 eye tracking dataset natively. In this dataset,
N = 30 (n = 15 during simultaneous fmri acquisition, n = 15
in a laboratory setting) participants watched the audiovisual movie Forrest Gump
while their eye movements were recorded with an Eyelink 1000. The movie was
presented in 8 segments of roughly 15 minutes of length. For all details on the
data acquisition, see the corresponding publication_ by Hanke and colleagues (2016).

The raw eyetracking data was classified into eye movements (fixations, saccades, 
postsaccadic oscillations, and pursuits) with the REMoDNaV (Robust Eye Movement
Detection for Natural Viewing) algorithm (Dar, Wagner & Hanke, in preperation).
These results can be found here_ and serve as input files for ``multimatch_forrest``.

Additionally, the studyforrest dataset contains extensive annotation. For
``multimatch_forrest``, the location-annotation_ (Häusler & Hanke, 2016) of the
movie is used to split the classified eye movement data into scanpaths of
user-specified length within a shot of the movie. This was implemented to take
into account that cuts in dynamic scenes generally lead to a strong center bias
in viewers (Carmi & Itti, 2006). The user can specify whether scanpaths should
start with the beginning of a shot, or, in order to include as little center
bias as possible, should be extraced to end precisely with the end of the shot
(and thus have the longest possible distance between shot onset and scanpath
onset).

The function will take two eyemovement datafiles of one run (one ~15 minute segment,
from two subjects respectively), annotation data of the corresponding run, and
an output path as required inputs and returns a .tsv event file. One row of the
file corresponds to one scanpath comparison, the columns are the onsets of the
compared scanpaths, the exact durations of the scanpaths, and similarity values
on the five dimensions per comparison.


.. _studyforrest: https://github.com/psychoinformatics-de/studyforrest-data-phase2
.. _here: https://github.com/psychoinformatics-de/studyforrest-data-eyemovementlabels
.. _publication: https://www.nature.com/articles/sdata201692
.. _location-annotation: https://github.com/psychoinformatics-de/studyforrest-data-annotations

execution via commandline
^^^^^^^^^^^^^^^^^^^^^^^^^

Just as multimatch_, ``multimatch_forrest`` also works easiest when executed
in a terminal as a single command line. The comparison of all scanpaths of the
default length (4.92 seconds, the median shot length of the movie) only needs the
command line keyword ``multimatch_forrest``, followed by two inputs, corresponding
to the remodnav_ outputs of two subjects in the same run,
the annotation file for the shots of the respecitve run and an output path.

.. code::

   multimatch_forrest path/to/sub-a_run-x.tsv path/to/sub-b_run-x.tsv
   path/to/shotannotation where/results/go

Additionally, the following options can be specified:

- ``--screensize``: in px, specify first x, than y dimensions. Default is 1280 x
  720px.
- ``--direction_threshold``: If two consecutive saccades have a small angle, they will be
  combined. Should be in degrees, such as ``45.0`` for 45°.
- ``--amplitude_threshold``: If two consecutive saccades are short, they will be
  combines. Should be in pixel, such as ``100.0`` for 100px.
- ``--duration_threshold``: Only if the intermediate fixation's durations are
  shorter than this threshold the above simplification will be performed. Should
  be in seconds, such as ``0.1`` for 100ms.

**Note**: If either direction- or amplitude threshold is specified as 0, no
grouping will be performed!

- ``--duration``: The approximate desired duration for a scanpaths in
  seconds, e.g. 3.5. Default: 4.92s (the median shotlength).
- ``--lduration``: Option to group short shots in the same locale (i.e no
  change in setting) together for longer scanpaths. Shots shorter than ``ldur``
  will be attempted to be grouped together.
- ``--position-offset``: Boolean, if True, scanpaths of ``dur`` length
  stop at shotoffset (instead of beginning at shot onset). Default: False.



.. _multimatch: https://multimatch.readthedocs.io/en/latest/multimatch.html
.. _remodnav: https://github.com/psychoinformatics-de/studyforrest-data-eyemovementlabels

execution within a python instance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you wish to use the functionality of multimatch within a running python
instance such as ipython, you can import the module and use the function
``docomparison``. Here is an example:

.. code::

   import numpy as np
   import pandas as pd
   import multimatch as m

   # read in necessary datafiles
   # annotations (run 1)
   shots = pd.read_csv('multimatch/tests/testdata/locations_run-1_events.tsv',
   sep='\t')

   # eyemovement data sub-10
   remodnav_1 = np.recfromcsv('multimatch/tests/testdata/sub-10_task-movie_run-1_events.tsv',
   delimiter='\t', dtype={'names': ('onset', 'duration', 'label', 'start_x',
   'start_y', 'end_x', 'end_y', 'amp', 'peak_vel', 'med_vel', 'avg_vel'),
   'formats': ('f8', 'f8', 'U10', 'f8', 'f8', 'f8','f8', 'f8', 'f8', 'f8', 'f8')})

   # eyemovement data sub-30
   remodnav_2 = np.recfromcsv('multimatch/tests/testdata/sub-30_task-movie_run-1_events.tsv',
   delimiter='\t', dtype={'names': ('onset', 'duration', 'label', 'start_x',
   'start_y', 'end_x', 'end_y', 'amp', 'peak_vel', 'med_vel', 'avg_vel'),
   'formats': ('f8', 'f8', 'U10', 'f8', 'f8', 'f8','f8', 'f8', 'f8', 'f8', 'f8')})

   # execute scanpath comparison
   similarities, onsets, durations = m.docomparison_forrest(shots, remodnav_1,
   remodnav_2, sz=[1280, 720], dur=3.0, ldur=0, offset=False,TDur=0, TAmp=0,
   TDir=0, grouping=False)




References
^^^^^^^^^^

Carmi, R. & Itti, L. (2006). Visual causes versus correlates of attentional
selection in dynamic scenes. *Vision Research*, 46, 4333 – 4345.

Hanke, M., Adelhöfer, N., Kottke, D., Iacovella, V., Sengupta, A., Kaule, F. R.,
Nigbur, R., Waite, A. Q., Baumgartner, F. & Stadler, J. (2016).
A studyforrest extension, simultaneous fMRI and eye gaze recordings during
prolonged natural stimulation. *Scientific Data*, 3:160092.

Häusler, C. O. & Hanke, M. (2016). An annotation of cuts, depicted locations,
and temproal progression in the motion picture “Forrest Gump”. *F1000Research*,
5:2273.

