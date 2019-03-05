**********************
An example compuation
**********************
The following section shows a multimatch use case to compute the scanpath
similarities of participants that watched the Hollywood movie Forrest Gump
during simultaneous fMRI acquisition.

Data and sample
^^^^^^^^^^^^^^^
Data for all analyses stems from the 2016 released extension of the studyforrest dataset
(Hanke et al., 2016; Sengupta et al., 2016). In this extension,
N = 15 right-handed participants (age range 21 - 39 years, mean age 29.4 years, six female,
normal or corrected-to-normal vision), who had previously participated in the studyforrest
project, watched the audio-visual movie 'Forrest Gump' (R. Zemeckis, Paramount Pictures, 1994)
during simultaneous fMRI and eye-tracking recording. The video track for the movie stimulus
was re-encoded from Blu-ray into H.264 video (1280 x 720px at 25 frames per second
(fps)). In accordance to the procedure in an earlier phase of the studyforrest project, the
movie was shortened by removing a few scenes less relevant for the major plot to keep
the fMRI recording session under two hours. The shortened movie was then split into
eight segments of roughly 15 minutes of length (for an overview on segment duration,
final stimulus content and detailed procedures see Hanke et al. (2014)).
Visual stimuli were projected on to a screen inside the bore of the magnet using
an LCD projector, and presented to the subjects through a front-reflective mirror on
top of the head coil at a viewing distance of 63cm. The screen dimensions were 26.5cm
x 21.2cm (corresponding to 1280 x 1024px) at a resolution of 720p at full width, with
a 60Hz video refresh rate (Sengupta et al., 2016). Eye-tracking was performed with
an Eyelink 1000 (software version 4.594) using monocular corneal reflection and pupil
tracking with a temporal resolution of eye gaze recordings of 1000Hz.
The camera was mounted at an approximate distance of 100cm to the left eye of subjects, which
was illuminated by an infrared light source (Hanke et al., 2016). Eye-tracking data were normalized such
that all gaze coordinates are in native movie frame pixels, with the top-left corner of
the movie frame located at (0, 0) and the lower-right corner located at (1280, 546)
(ibid.). The amount of unusable data, primarily due to signal loss during eye blinks,
ranged from less than 1 to 15% for 13 of the 15 in-scanner subjects (the other two
subjects' data contained 85 and 36% of data loss, respectively). In-scanner acquisition
had an approximate spatial uncertainty of 40px according to the calibration procedure
(ibid.).

.. _studyforrest: https://github.com/psychoinformatics-de/studyforrest-data-phase2

Event detection and scanpath derivation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Raw gaze data was classified into different categories of eye movements
with an adaptive, data-driven algorithm for robust eye movement detection for natural
viewing (REMoDNaV_ ) in Python. The algorithm categorizes the raw data into
saccades, fixations, smooth pursuits, and post-saccadic oscillations
(glissades), and disregards any unclassifiable data (such as blinks). It was specifically
developed to compute robust results even under high noise conditions.
For an overview of the algorithmic details and evaluation of REMoDNaV compared to
contemporary algorithms and human annotations, please see the respective publication_
(Dar et al., in preparation) or take a look at the REMoDNaV_ module.
The eye events are reported together with their start and end coordinates, their onsets
and durations in seconds, their velocity, and the average pupil size.
Fixation vectors as input for multimatch were derived from the REMoDNaV output.
As the stimulus was dynamic with moving targets that evoke smooth pursuit movements,
such pusuit events are categorized to be
an eye movement category of 'visual intake', just as fixation. Therefore, in a first step,
the start and end of pursuit movements were included in scanpaths to compare as well.
In a second step, the continuous eye movement data (~15 min per run) was split into shots
corresponding to segments that did not contain scene changes between depicted
locales using the published location annotation for the movie (Häusler & Hanke,
2016). This was done to accommodate the fact that subjects gazes have
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
In total, 533 scanpaths were extracted from the movie. The median duration of extracted scanpath
duration was 4.39 seconds (mean = 4.36s).
The following figures give an overview of the similarity computations.
Figures 1 and 2 display a frame within the segments in the first run of the movie
with the lowest and highest group-level similarity (averaged across the five dimensions).

 .. figure:: ../img/low_sim.png
   :scale: 50%
   :alt: low similarity segment

   One frame from the segment within the first run of the movie with the **lowest** average group-level similarity.
   The circles represent participants center of eye gaze.

 .. figure:: ../img/max_sim.png
   :scale: 50%
   :alt: high similarity segment

   One frame from the segment within the first run of the movie with the **highest** average group-level similarity.
   The circles represent participants center of eye gaze.


The overall similarity of gaze was high, however, there were consistent differences between
dimensions. The Shape, Length and Position
dimension displayed very high similarities, and the average Duration similarity was
the lowest of all dimensions.
Medians and means correspond closely, and standard
deviations are very small. This is also highlighted by Figure 3.

   =========   ===========  =========
   Variable    mean [SD]    median
   =========   ===========  =========
   Shape       0.97 [0.01]  0.97
   Position    0.88 [0.03]  0.89
   Length      0.96 [0.01]  0.96
   Duration    0.54 [0.05]  0.55
   Direction   0.72 [0.05]  0.71
   =========   ===========  =========


 .. figure:: ../img/sim_per_dimension.png
   :scale: 100%
   :alt: distribution of similarity measures

   Distribution of similarity measures throughout the movie. Note the extremely high
   position and length dimension.

Discussion
^^^^^^^^^^

As evident from the previous table and figure, scanpaths were almost
perfectly similar on the dimensions vector length and vector position.
This is likely at least partially due to the scanpath alignment based on the scanpath shape.
Scanpaths were also highly similar on the position dimension, which demonstrates a strong
gaze control of the movie stimulus. Subjects scanpaths differed more substantially on
the dimensions direction and duration, which indicates differences in fixation dwelling times and saccadic angle. Thus, the general points of interest (as evident from high
similarities in position, length and shape) were similar across subject, but dierences in
direction and duration might indicate interindividually dierent exploration strategies.
All dimensions show a remarkable consistency in similarity measures as evident from
the small standard deviations. This might indicate a consistently high level of exogenous
attentional control by the movie stimulus. This finding is consistent with research on
viewing behavior during movies: Unlike during static image viewing, the spatio-temporal
gaze behavior of multiple viewers exhibits a substantial degree of coordination in movie
watching. Smith and Henderson (2008) cued the term *attentional synchrony* for this
phenomenon. During attentional synchrony, viewers gazes cluster around a small portion
of the screen at any one moment. Goldstein et al. (2007), for example, found the
distribution of fixations of viewers to occupy less than 12% of the total screen area
in more than 50% of the time in six Hollywood movies. In a comparison between
different types of static and dynamic visual stimuli, Dorr et al. (2010) found the
highest consistency between viewers eyegazes during professionally produced (Hollywood)
movies, likely largely due to the use of cinematic composition of scenes, deliberate
camera work and editing. Hasson et al. (2008) found high correspondence in gaze behavior
across subjects, even for backwards presentations of movies.

The results obtained with the multimatch algorithm from the Hollywood movie
Forrest Gump, therefore, are consistent with known properties of gaze behavior
during movie watching. This analysis has furthermore demonstrated one way of using
multimatchs scanpath comparison on a grouplevel similarity computation per segment.
If you have any questions about this example, please ask here_.

 .. _here: https://github.com/AdinaWagner/multimatch/issues/new



