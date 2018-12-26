[![Build Status](https://travis-ci.com/AdinaWagner/multimatch.svg?branch=master)](https://travis-ci.com/AdinaWagner/multimatch)
[![codecov](https://codecov.io/gh/AdinaWagner/multimatch/branch/master/graph/badge.svg)](https://codecov.io/gh/AdinaWagner/multimatch)
[![Documentation](https://readthedocs.org/projects/multimatch/badge/?version=latest)](https://multimatch.readthedocs.io/en/latest/)
[![PyPIversion](https://badge.fury.io/py/multimatch.svg)](https://badge.fury.io/py/multimatch)
 [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
# multimatch
## Reimplementation of MultiMatch toolbox (Dewhurst et al., 2012) in Python.

The **MultiMatch** method proposed by Jarodzka, Holmqvist and Nyström (2010),
implemented in Matlab as the MultiMatch toolbox and validated by Dewhurst
and colleagues (2012) is a vector-based, multi-dimensional approach to
compute scanpath similarity.

The method represents scanpaths as geometrical vectors in a two-dimensional
space: Any scanpath is build up of a vector sequence in which the vectors
represent saccades, and the start and end position of saccade vectors represent
fixations. Two such sequences (which can differ in length) are compared on the
five dimensions **'vector shape'**, **'vector length'** (saccadic amplitude),
**'vector position'**, **'vector direction'** and **'fixation duration'** for a
multidimensional similarity evaluation. The original Matlab toolbox was kindly
provided via email by Dr. Richard Dewhurst and the method was ported into Python
with the intent of providing an open source alternative to the matlab toolbox.
Additionally, it contains options to compute scanpath similarities from the
studyforrest phase 2 eyetracking dataset, in which participants (n = 15 during
fmri acquisition, n = 15 in lab) watched the movie 'Forrest Gump'.

### Installation instructions

It is recommended to use a dedicated virtualenv:

    # create and enter a new virtual environment (optional)
    virtualenv --python=python3 ~/env/multimatch
    . ~/env/multimatch/bin/activate

multimatch can be installed via pip. To automatically install multimatch with all
dependencies, use:

    # install from pyPi
    pip install multimatch


### Method overview

The method takes two n x 3 fixation vectors (x-coordinate, y-coordinate,
duration) of two scanpaths as its input.

- **Step 1: Representation of scanpaths as vector sequences**
An idealized saccade is represented as the shortest distance between two
fixations. The Cartesian coordinates of the fixations are thus the starting and
ending points of a saccade. The length of a saccade in x direction is computed
as the difference in x coordinates of starting and ending point. The length of a
saccade in y direction is computed accordingly. To represent a saccade as a
vector in two-dimensional space, the lengths in x and y directions are
transformed into polar coordinates (length from coordinate origin (Rho), polar
angle in radians (Theta)) by means of trigonometry.

- **Step 2: Scanpath simplification**
Scanpaths are simplified based on angle and amplitude (length) to reduce their
complexity. Two or more saccades are grouped together if angles between two
consecutive saccades are below an angular threshold ```TAmp```, and intermediate
fixations are shorter than a duration threshold ```TDur```, or if the amplitude
of successive saccades is below a length threshold ```TAmp``` and the
surrounding fixation duration. As such, small, locally contained saccades, and
saccades in the same general direction are summed to form larger, less complex
saccades (Dewhurst et al., 2012). This process is repeated
until no further simplifications are made.
Thresholds can be set according to use case. The original simplification algorithm
implements an angular threshold of 45° and an amplitude threshold of 10% of the
screen diagonal (Jarodzka, Holmqvist & Nyström, 2010).

- **Step 3: Temporal Alignment**
Two simplified scanpaths are temporally aligned in order to find pairings of
saccade vectors to compare. The aim is not necessarily to align two saccade
vectors that constitute the same component in  their respective vector sequence,
but those two vectors that are the most similar while preserving temporal order.
In this way, a stray saccade in one of the two scanpaths does not lead to an
overall low similarity rating, and it is further possible to compare scanpaths
of unequal length.  To do so, all possible pairings of saccades are evaluated in
similarity by their shape (i.e. vector differences). More formally, the vector
difference between each element i in scanpath S1 = {u1, u2, …, um} and each
element j in scanpath S2 = {v1, v2, …, vn} is computed and stored in Matrix M
as a weight. Low weights correspond to high similarity. An adjacency matrix of
size M is build, defining rules on which connection between matrix elements are
allowed: In order to take temporal sequence of saccades into account, connections
can only be made to the right, below or below-right.
Together, matrices M and the adjacency matrix constitute a matrix representation
of a directed, weighted graph. The elements of the matrix are the
nodes, the connection rules constitute edges and the weights define the cost
associated with each connection.

- **Step 4: Scanpath Selection**
A Dijkstra algorithm (Dijksta, 1959) is used to find the shortest path from the
the first two saccade vectors to the last two saccade vectors. “Shortest” path
is defined as the connection between nodes with the lowest possible sum of
weights.

- **Step 5: Similarity Calculation**
Five measures of scanpath similarity are computed on the aligned scanpaths. This
is done by performing simple vector arithmetic on all aligned saccade pairs
(u_i, v_j), taking the median of the results and normalizing it. As a result,
all five measures are in range [0, 1] with higher values indicating higher
similarity between scanpaths on the given dimension (Anderson, Anderson,
Kingstone & Bischof, 2015).

For more details on the original algorithm, please see Dewhurst et al. (2012).

### Examplary usage of multimatch in a terminal

**required inputs:**
- two tab-separated files with nx3 fixation vectors (x, y, duration)

`` multimatch data/fixvectors/segment_10_sub-19.tsv data/fixvectors/segment_10_sub-01.tsv ``



**optional inputs:**
- --screensize: in pixel, supply first x and then y dimension. The default size is 1280 x 720px

`` multimatch data/fixvectors/segment_10_sub-19.tsv data/fixvectors/segment_10_sub-01.tsv --screensize 1280 720 ``

if scanpath simplification should be performed, please specify in addition
- --amplitude_threshold (-am) in px
- --direction_threshold (-di) in degree
- --duration_threshold (-du) in seconds

Example usage with grouping:

`` multimatch data/fixvectors/segment_10_sub-19.tsv
data/fixvectors/segment_10_sub-01.tsv --direction_threshold 45.0
--duration_threshold 0.3 --amplitude_threshold 147.0 ``


## multimatch_forrest

The package further contains an implementation of the method specifically for
use with the studyforrest phase 2 eyetracking dataset. Inputs should be the
eyetracking data classified into eye-movement events by remodnav
(https://github.com/psychoinformatics-de/remodnav).
As such, inputs are two eye-movement classified tsv-files of a given
run of studyforrest and the respective location annotation
(https://github.com/psychoinformatics-de/studyforrest-data-annotations) file as
well an output path. It will create scanpaths from the shots of the movie stimlus
and compute the necessary nx3 data structure on its own.

### Examplary usage of multimatch_forrest in a terminal:

``multimatch_forrest
multimatch/tests/testdata/sub-10_task-movie_run-1_events.tsv multimatch/tests/testdata/sub-30_task-movie_run-1_events.tsv multimatch/tests/testdata/studyforrest-data-annotations/segments/avmovie/locations_run-1_events.tsv output/run-1/sub-10vssub-30``


### References:

Dewhurst, R., Nyström, M., Jarodzka, H., Foulsham, T., Johansson, R. &
Holmqvist, K. (2012). It depends on how you look at it: scanpath comparison in
multiple dimensions with MultiMatch, a vector-based approach. Behaviour Research
Methods, 44(4), 1079-1100.

Dijkstra, E. W. (1959). A note on two problems in connexion withgraphs.
Numerische Mathematik, 1, 269–271.

Jarodzka, H., Holmqvist, K., & Nyström, M. (eds.) (2010). A vector-based,
multidimensional scanpath similarity measure. In Proceedings of the 2010
symposium on eye-tracking research & applications (pp. 211-218). ACM.

