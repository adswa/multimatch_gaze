[![Build Status](https://travis-ci.com/adswa/multimatch_gaze.svg?branch=master)](https://travis-ci.com/adswa/multimatch_gaze)
[![codecov](https://codecov.io/gh/adswa/multimatch_gaze/branch/master/graph/badge.svg)](https://codecov.io/gh/adswa/multimatch_gaze)
[![Documentation](https://readthedocs.org/projects/multimatch/badge/?version=latest)](https://multimatch.readthedocs.io/en/latest/)
[![PyPI version](https://badge.fury.io/py/multimatch-gaze.svg)](https://badge.fury.io/py/multimatch-gaze)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build status](https://ci.appveyor.com/api/projects/status/wrphckxqjrfut703?svg=true)](https://ci.appveyor.com/project/adswa/multimatch_gaze)
[![DOI](https://zenodo.org/badge/151181532.svg)](https://zenodo.org/badge/latestdoi/151181532)

# multimatch-gaze
## Reimplementation of MultiMatch toolbox (Dewhurst et al., 2012) in Python.

The **MultiMatch** method proposed by Jarodzka, Holmqvist and Nyström (2010),
implemented in Matlab as the MultiMatch toolbox and validated by Dewhurst
and colleagues (2012) is a vector-based, multi-dimensional approach to
compute scan path similarity.

For a complete overview of this software, please take a look at the
[Documentation](https://multimatch.readthedocs.io/en/latest)

The method represents scan paths as geometrical vectors in a two-dimensional
space: Any scan path is build up of a vector sequence in which the vectors
represent saccades, and the start and end position of saccade vectors represent
fixations. Two such sequences (which can differ in length) are compared on the
five dimensions **'vector shape'**, **'vector length'** (saccadic amplitude),
**'vector position'**, **'vector direction'** and **'fixation duration'** for a
multidimensional similarity evaluation (all in range [0, 1] with 0 denoting
maximal dissimilarity and 1 denoting identical scan paths on the given measure).
The original Matlab toolbox was kindly
provided via email by Dr. Richard Dewhurst and the method was ported into Python
with the intent of providing an open source alternative to the matlab toolbox.

### Installation instructions

It is recommended to use a dedicated virtualenv:

    # create and enter a new virtual environment (optional)
    virtualenv --python=python3 ~/env/multimatch
    . ~/env/multimatch/bin/activate

multimatch-gaze can be installed via pip. To automatically install multimatch-gaze with all
dependencies, use:

    # install from pyPi
    pip install multimatch-gaze


### Support/Contributing

Bug reports, feedback, or any other contribution are always appreciated. To
report a bug, request a feature, or ask a question, please open an
[issue](https://github.com/adswa/multimatch_gaze/issues/new).
[Pull requests](https://help.github.com/en/articles/creating-a-pull-request-from-a-fork)
are always welcome.


### Examplary usage of multimatch-gaze in a terminal

**required inputs:**
- two tab-separated files with nx3 fixation vectors (x coordinate in px, y coordinate in px, duration)
- screensize in px (x dimension, y dimension)

`` multimatch-gaze data/fixvectors/segment_10_sub-19.tsv data/fixvectors/segment_10_sub-01.tsv 1280 720 ``



**optional inputs:**

if scan path simplification should be performed, please specify in addition
- --amplitude-threshold (-am) in px
- --direction-threshold (-di) in degree
- --duration-threshold (-du) in seconds

Example usage with grouping:

`` multimatch-gaze data/fixvectors/segment_10_sub-19.tsv
data/fixvectors/segment_10_sub-01.tsv 1280 720 --direction-threshold 45.0
--duration-threshold 0.3 --amplitude-threshold 147.0 ``

**REMoDNaV helper:**

Eye movement event detection results produced by [REMoDNaV](https://github.com/psychoinformatics-de/remodnav)
can be read in natively by multimatch-gaze. To indicate that datafiles are REMoDNaV outputs, supply the
``--remodnav`` parameter.

`` multimatch-gaze data/remodnav_samples/sub-01_task-movie_run-1_events.tsv
data/remodnav_samples/sub-01_task-movie_run-2_events.tsv 1280 720 --remodnav ``

REMoDNaV can classify smooth pursuit movements. As a consequence, when using REMoDNaV output, users need to
indicate how these events should be treated. By default, multimatch-gaze will discard pursuits. In some
circumstances, however, it can be useful to include pursuit information. Moving stimuli for example would
evoke a pursuit movement during visual intake. When specifying the ``--pursuit keep`` parameter, the start
and end points of pursuits will be included in the scan path.

`` multimatch-gaze data/remodnav_samples/sub-01_task-movie_run-1_events.tsv
data/remodnav_samples/sub-01_task-movie_run-2_events.tsv 1280 720 --remodnav --pursuit keep``


### References:

Dewhurst, R., Nyström, M., Jarodzka, H., Foulsham, T., Johansson, R. &
Holmqvist, K. (2012). It depends on how you look at it: scanpath comparison in
multiple dimensions with MultiMatch, a vector-based approach. Behaviour Research
Methods, 44(4), 1079-1100. [doi: 10.3758/s13428-012-0212-2.](https://doi.org/10.3758/s13428-012-0212-2)

Dijkstra, E. W. (1959). A note on two problems in connexion withgraphs.
Numerische Mathematik, 1, 269–271. [https://doi.org/10.1007/BF01386390](https://doi.org/10.1007/BF01386390)

Jarodzka, H., Holmqvist, K., & Nyström, M. (eds.) (2010). A vector-based,
multidimensional scanpath similarity measure. In Proceedings of the 2010
symposium on eye-tracking research & applications (pp. 211-218). ACM.
[doi: 10.1145/1743666.1743718](https://doi.org/10.1145/1743666.1743718)

