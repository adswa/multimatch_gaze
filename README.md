# MultiMatch
Reimplementation of MultiMatch toolbox (Dewhurst et al., 2012) in Python.

The MultiMatch method proposed by Jarodzka, Holmqvist and Nyström (2010),
implemented in Matlab as the MultiMatch toolbox and validated by Dewhurst
and colleagues (2012) is a vector-based, multi-dimensional approach to
measuring scanpath similarity.

The method represents scanpaths as geometrical vectors in a two-dimensional
space: Any scanpath is build up of a vector sequence in which the vectors
represent saccades, and the start and end position of saccade vectors represent
fixations. Two such sequences (which can differ in length) are compared on the
five dimensions 'vector shape', 'vector length' (saccadic amplitude), 'vector
position', 'vector direction' and 'fixation duration' for a multidimensional
similarity evaluation. The original Matlab toolbox was kindly provided via email
by Dr. Richard Dewhurst and the method was ported into Python with the intent of
comparing scanpaths evoked by movie 'Forrest Gump' in subjects.

The orginal method takes n x 3 fixation vectors of two scanpaths with x- and
y-coordinates of the fixations and their duration as variables as its input. In
a first step, based on the coordinates and durations of fixations, the scanpaths
are represented as a vector sequences:
An idealized saccade is represented as the shortest distance between two
fixations. The Cartesian coordinates of the fixations are thus the starting and
ending points of a saccade. The length of a saccade in x direction is computed
as the difference in x coordinates of starting and ending point. The length of a
saccade in y direction is computed in the same way but with y coordinates,
respectively. To represent a saccade as a vector in two-dimensional space, the
lengths in x and y directions are transformed into polar coordinates, namely
length from coordinate origin (Rho) and polar angle in radians (Theta) by means
of trigonometry.

In a second step, the scanpaths are simplified according to angle and amplitude
(length) to reduce their complexity. Two or more saccades are grouped together
if angles between two consecutive saccades are below an angular threshold, or if
the amplitude of successive saccades is below a length threshold. As such,
small, locally contained saccades, and saccades in the same general direction
are summed to form larger, less complex saccades (Dewhurst et al., 2012).
Thresholds can be set according to use case, but the simplification algorithm
implements an angular threshold of 45° and an amplitude threshold of 10% of the
screen diagonal (Jarodzka, Holmqvist & Nyström, 2010). This process is repeated
until no further simplifications are made.

In a third step, the two simplified scanpaths are temporally aligned in order to
find pairings of saccade vectors to compare. The aim is to not necessarily align
two saccade vectors that constitute the same component in  their respective
vector sequence, but those two vectors that are the most similar while still
preserving temporal order. In this way, a stray saccade in one of the two
saccades does not lead to an overall low similarity rating, and it is further
possible to compare scanpaths of unequal length.  To do so, all possible
pairings of saccades are evaluated in similarity by a metric such as shape (i.e.
vector differences between all pairings). More formally, the vector difference
between each element i in scanpath S1 = {u1, u2, …, um} and each element j in
scanpath S2 = {v1, v2, …, vn} is computed an stored in a Matrix M. The
components of matrix M are the weights that denote similarity. Low weights
correspond to high similarity. In a next step, an adjacency matrix of size M is
build. The adjacency matrix defines rules on which connection between matrix
elements are allowed: In order to take temporal sequence of saccades into
account, connections can only be made to the right, below or below-right.
Together, matrices M and the adjacency matrix constitute a matrix representation
of a directed, weighted graph. The elements of the matrix are the
nodes, the connection rules constitute edges and the weights define the cost
associated with each connection.

A Dijkstra algorithm (Dijksta, 1959) is used to find the shortest path from the
top left node, the first two saccade vectors, to the bottom right node, the last
two saccade vectors. “Shortest” path is defined as the connection between nodes
with the lowest possible sum of weights. This algorithm is common in way-finding
and can be compared to a simple car navigation system that finds the shortest
driving distance between two cities. The path returned by the Dijkstra algorithm
is a sequence of indexes, denoting pairings of saccade vectors from each
scanpath, and as such the desired alignment of scanpaths (Dewhurst et al.,
2012).

Finally, in a last step, five measures of scanpath similarity are computed. This
is done by performing simple vector arithmetic on all aligned saccade pairs
(u_i, v_j), taking the average of the results and normalizing this average
according to a certain metric. As a result, all five measures are in range [0,
1] with higher values indicating higher similarity between scanpaths on the
given dimension (Anderson, Anderson, Kingstone & Bischof, 2015).
For vector shape similarity, the vector difference is computed as  u_i – v_j
and normalized by 2x the screen diagonal (the maximum theoretic value).
For vector length similarity, the difference in length between aligned saccade
vectors is computed and normalized by the screen diagonal.
For vector direction similarity, the angular distance in radians between aligned
saccade vectors is computed and normalized by pi.
For fixation position similarity, the difference in fixation coordinates is
calculated as the Euclidean distance and normalized by the screen diagonal.
For fixation duration similarity, the difference in fixation durations between
aligned fixations is computed and normalized against the maximum duration of the
two fixation durations being compared (Dewhurst et al., 2012).

For more details on the original algorithm, please see Dewhurst et al. (2012).

As of now, the Python code in this repository is custom made for use with files
stemming from the studyforrest dataset's eyetracking data classified into
eyemovement events by remodnav (https://github.com/psychoinformatics-de/remodnav).
As such, the script will take two eye-movement classified tsv-files of a given
run of studyforrest and the respective location annotation
(https://github.com/psychoinformatics-de/studyforrest-data-annotations) file as
well an output path as inputs. It will create scanpaths from the shots of the
movie stimlus and compute the necessary nx3 data structure on its own.
In the future, the code will hopefully be updated to be a standalone used with
ready-to-use nx3 fixation vectors with the additional possibility of usage
within the studyforrest context.

Examplary usage of the script in a terminal:

./MultiMatch.py -i inputs/eyemovements/sub-01/sub-01_task-movie_run-1_events.tsv
-j inputs/eyemovements/sub-02/sub-02_task-movie_run-1_events.tsv -k
inputs/studyforrest-data-annotations/segments/avmovie/locations_run-1_events.tsv
-o output/run-1/sub-01vssub-02




References:

Dewhurst, R., Nyström, M., Jarodzka, H., Foulsham, T., Johansson, R. &
Holmqvist, K. (2012). It depends on how you look at it: scanpath comparison in
multiple dimensions with MultiMatch, a vector-based approach. Behaviour Research
Methods, 44(4), 1079-1100.

Dijkstra, E. W. (1959). A note on two problems in connexion withgraphs.
Numerische Mathematik, 1, 269–271.

Jarodzka, H., Holmqvist, K., & Nyström, M. (eds.) (2010). A vector-based,
multidimensional scanpath similarity measure. In Proceedings of the 2010
symposium on eye-tracking research & applications (pp. 211-218). ACM.

