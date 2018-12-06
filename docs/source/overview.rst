*********
overview
*********

multimatch is a python-based reimplementation of Dewhursts and colleagues (2012)
implementation of the MultiMatch algorithm (Jarodzka, Holmqvist & Nyström, 2010)
in a matlab toolbox. It is a vector-based, multi-dimensional approach to
compute scanpath similarity.

The method represents scanpaths as geometrical vectors in a two-dimensional
space: Any scanpath is build up of a vector sequence in which the vectors
represent saccades, and the start and end position of saccade vectors represent
fixations. Two such sequences (which can differ in length) are compared on the
five dimensions **vector shape**, **vector length** (saccadic amplitude), **vector
position**, **vector direction** and **fixation duration** for a multidimensional
similarity evaluation.

 .. figure:: ../img/dimensions.png
   :scale: 100%
   :alt: scanpath dimensions

   Dimensions of scanpath comparison, taken from Dewhurst et al., 2012

The original Matlab toolbox was kindly provided via email by Dr. Richard Dewhurst
and the method was ported into Python with the intent of providing an open source
alternative to the matlab toolbox.

The module provides the possibility to compute the similarity of two scanpaths
with a terminal command or within a python instance (see section multimatch_).

 .. _multimatch: https://multimatch.readthedocs.io/en/latest/multimatch.html

Additionally, the reimplementation contains
options to compute scanpath similarities from the studyforrest phase 2
eyetracking dataset, in which participants (n = 15 during fmri acquisition, n =
15 in lab) watched the movie 'Forrest Gump'. These additional functions of the
multimatch module are mainly concerned with automatic scanpath extraction from
the approximately 15 minutes spanning and several eye-movement categories
containing eye-tracking event files by user-defined rules. Inputs for the
studyforrest-version stem from remodnav1 (Dar, Wagner & Hanke, in preparation),
an adaptive eye-event classification algorithm used on the studyforrest
eye-tracking data and publicly available from
github.com/psychoinformatics-de/studyforrest-data-eyemovementlabels. For details
on how to use this functionality, see section multimatch_forrest_.

 .. _multimatch_forrest: https://multimatch.readthedocs.io/en/latest/multimatch_forrest.html


References
^^^^^^^^^^
Dewhurst, R., Nyström, M., Jarodzka, H., Foulsham, T., Johansson, R. &
Holmqvist, K. (2012). It depends on how you look at it: scanpath comparison in
multiple dimensions with MultiMatch, a vector-based approach. Behaviour Research
Methods, 44(4), 1079-1100.

Jarodzka, H., Holmqvist, K., & Nyström, M. (eds.) (2010). A vector-based,
multidimensional scanpath similarity measure. In Proceedings of the 2010
symposium on eye-tracking research & applications (pp. 211-218). ACM.
