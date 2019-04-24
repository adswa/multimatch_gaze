[![Build Status](https://travis-ci.com/adswa/multimatch.svg?branch=master)](https://travis-ci.com/adswa/multimatch)
[![codecov](https://codecov.io/gh/adswa/multimatch/branch/master/graph/badge.svg)](https://codecov.io/gh/adswa/multimatch)
[![Documentation](https://readthedocs.org/projects/multimatch/badge/?version=latest)](https://multimatch.readthedocs.io/en/latest/)
[![PyPIversion](https://badge.fury.io/py/multimatch.svg)](https://badge.fury.io/py/multimatch)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build status](https://ci.appveyor.com/api/projects/status/wrphckxqjrfut703?svg=true)](https://ci.appveyor.com/project/adswa/multimatch)


# multimatch
## Reimplementation of MultiMatch toolbox (Dewhurst et al., 2012) in Python.

The **MultiMatch** method proposed by Jarodzka, Holmqvist and Nyström (2010),
implemented in Matlab as the MultiMatch toolbox and validated by Dewhurst
and colleagues (2012) is a vector-based, multi-dimensional approach to
compute scanpath similarity.

For a complete overview of this software, please take a look at the
[Documentation](https://multimatch.readthedocs.io/en/latest)

The method represents scanpaths as geometrical vectors in a two-dimensional
space: Any scanpath is build up of a vector sequence in which the vectors
represent saccades, and the start and end position of saccade vectors represent
fixations. Two such sequences (which can differ in length) are compared on the
five dimensions **'vector shape'**, **'vector length'** (saccadic amplitude),
**'vector position'**, **'vector direction'** and **'fixation duration'** for a
multidimensional similarity evaluation (all in range [0, 1] with 0 denoting
maximal dissimilarity and 1 denoting identical scanpaths on the given measure).
The original Matlab toolbox was kindly
provided via email by Dr. Richard Dewhurst and the method was ported into Python
with the intent of providing an open source alternative to the matlab toolbox.

### Installation instructions

It is recommended to use a dedicated virtualenv:

    # create and enter a new virtual environment (optional)
    virtualenv --python=python3 ~/env/multimatch
    . ~/env/multimatch/bin/activate

multimatch can be installed via pip. To automatically install multimatch with all
dependencies, use:

    # install from pyPi
    pip install multimatch


## Support/Contributing

Bug reports, feedback, or any other contribution are always appreciated. To
report a bug, request a feature, or ask a question, please open an
[issue](https://github.com/adswa/multimatch/issues/new).
[Pull requests](https://help.github.com/en/articles/creating-a-pull-request-from-a-fork)
are always welcome.


### Examplary usage of multimatch in a terminal

**required inputs:**
- two tab-separated files with nx3 fixation vectors (x coordinate in px, y coordinate in px, duration)

`` multimatch data/fixvectors/segment_10_sub-19.tsv data/fixvectors/segment_10_sub-01.tsv ``



**optional inputs:**
- --screensize: in pixel, supply first x and then y dimension. The default size is 1280 x 720px

`` multimatch data/fixvectors/segment_10_sub-19.tsv data/fixvectors/segment_10_sub-01.tsv --screensize 1280 720 ``

if scanpath simplification should be performed, please specify in addition
- --amplitude-threshold (-am) in px
- --direction-threshold (-di) in degree
- --duration-threshold (-du) in seconds

Example usage with grouping:

`` multimatch data/fixvectors/segment_10_sub-19.tsv
data/fixvectors/segment_10_sub-01.tsv --direction-threshold 45.0
--duration-threshold 0.3 --amplitude-threshold 147.0 ``


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

