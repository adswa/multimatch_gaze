*********
overview
*********

**multimatch** is a python-based reimplementation of Dewhursts and colleagues (2012)
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

For a more detailed overview of the algorithm, take a look at the method_
section, or take a look at the original publication by Dewhurst_ et al. (2012)
and Jarodzka_ et al. (2010).

.. _method: https://multimatch.readthedocs.io/en/latest/method.html

.. _Dewhurst: https://link.springer.com/article/10.3758%2Fs13428-012-0212-2

.. _Jarodzka: http://portal.research.lu.se/ws/files/5608175/1539210.PDF

The original Matlab toolbox was kindly provided via email by Dr. Richard Dewhurst
and the method was ported into Python with the intent of providing an open source
alternative to the matlab toolbox.


Functionality
^^^^^^^^^^^^^

The module provides the possibility to compute the similarity of two scanpaths
with a terminal command or within a python instance (see section multimatch_).

 .. _multimatch: https://multimatch.readthedocs.io/en/latest/multimatch.html


References
^^^^^^^^^^
Dewhurst, R., Nyström, M., Jarodzka, H., Foulsham, T., Johansson, R. &
Holmqvist, K. (2012). It depends on how you look at it: scanpath comparison in
multiple dimensions with MultiMatch, a vector-based approach. *Behaviour Research
Methods, 44* (4), 1079-1100.

Jarodzka, H., Holmqvist, K., & Nyström, M. (eds.) (2010). A vector-based,
multidimensional scanpath similarity measure. In *Proceedings of the 2010
symposium on eye-tracking research & applications* (pp. 211-218). ACM.
