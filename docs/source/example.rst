**********************
An example compuation
**********************
The following section shows a multimatch use case to compute the scanpath
similarities of participants that watched the Hollywood movie Forrest Gump
during simultaneous fMRI acquisition.

Data and sample
^^^^^^^^^^^^^^^
Data stems from the 2016 release of the studyforrest_ dataset. N = 15 subjects
watched the audiovisual movie "Forrest Gump". The dataset contains participants
eye gaze location (recorded with an Eyelink 1000).

.. _studyforrest: https://github.com/psychoinformatics-de/studyforrest-data-phase2

Preprocessing
^^^^^^^^^^^^^
Raw eye gaze data were normalized such that all gaze coordinates are in native
movie frame pixels, with the top-lft corner of the movie frame located at (0,0)
and the lower right corner located at (1280, 546) (Hanke et al., 2016).
Subsequently, the raw gaze data was classified into different categories of eye movements
with a data-driven algorithm for robust eye movement detection for natural
viewing (REMoDNaV_ ) in Python. The algorithm categorizes the raw data into
saccades, fixations, smooth pursuits, and post-saccadic oscillations
(glissades), and disregards any unclassifiable data (such as blinks). The eye
events are reported together with their start and end coordinates, their onsets
and durations in seconds, their velocity, and the average pupil size. For
further information regarding the algorithm and its use, please see its
corresponding publication_.
Fixation vectors were derived from the REMoDNaV output. As the stimulus was
dynamic, in a first step, the start and end of pursuit movements were relabled as fixations as
well. This was done to accommodate the fact that most areas of interest in a
movie are in motion.
In a second step, the continuous eye movement data was split into shots
corresponding to segments that did not contain scene changes between depicted
locales using the published location annotation for the movie (Häusler & Hanke,
2016). This was done to accomadate the fact that subjects gazes have
a bias towards the center in Hollywood movies (Tseng et al. 2009). This bias can
at least in part be traced back to a strong center bias directly after cuts in
dynamic scenes. Lastly, within each segment, scanpaths of the median shot length
of ~4.92 seconds. To further evade any problems associated with the center bias,
scanpaths were extracted from the end of the segment: The last oculomotor event
within the range of the segment marked the end of a scanpath. As such, scanpaths
began maximally distant to the snippet onset.


.. _REmoDNaV: https://github.com/psychoinformatics-de/remodnav
.. _publication: pathtoourpublication.de

multimatch application
^^^^^^^^^^^^^^^^^^^^^^
Overall scanpath similarities were computed in a two-step procedure. First,
scanpath comparisons of all scanpaths from the same shot of two subjects were
calculated for all possible pairs of subject. This resulted in 105 combinations
for N = 15 subjects. These comparisons were done without any further
simplification (i.e. no use of the direction, length, and duration thresholds),
as even minor differences in scanpaths obtained from a movie can correspond to
major differences in attended visual stimuli. In a second step, the resulting
similarities for each of the five similarity dimensions were averaged. Thus, for
each snippet longer than 4.92s five similarity measures were computed that
represented the average similarity of scanpaths of all subjects on the given
dimension.
The results of this computation can be found on Github_.

.. _Github: https://www.github.com/AdinaWagner/multimatch_forrest

Results
^^^^^^^

The following graphics give an overview of the similarity computations.


