***************
multimatch-gaze
***************


**multimatch-gaze** is a python-based reimplementation of the MultiMatch algorithm
(Jarodzka, Holmqvist & Nyström, 2010).
The original Matlab toolbox was kindly provided via email by Dr. Richard Dewhurst
and the method was ported into Python with the intent of providing an open source
alternative to the matlab toolbox.


The module provides the possibility to compute the similarity of two scanpaths
with a terminal command or within a python instance (see section API_).

 .. _API: https://multimatch_gaze.readthedocs.io/en/latest/api.html



==============
gettingstarted
==============

Installation
************

It is recommended to use a dedicated virtualenv_.

.. _virtualenv: https://virtualenv.pypa.io

.. code::

   # create and enter a new virtual environment (optional)
   virtualenv --python=python3 ~/env/multimatch
   . ~/env/multimatch/bin/activate


Via pip install
---------------


multimatch-gaze can be installed via pip_ (**P**\ip **I**\nstalls **P**\ython). To
automatically install multimatch-gaze with all dependencies type::

   pip install multimatch-gaze

.. _pip: https://pip.pypa.io


Via Github
----------

The source code for multimatch-gaze can be found on Github_.

.. _Github: https://github.com/adswa/multimatch_gaze


A short tutorial...
-------------------

... session to get a first hang on how to use multimatch-gaze can be
conducted by cloning the Github repository and executing the
examples provided in the API_
section. The data used in these examples corresponds to the
data found in the repository.

.. _API: https://multimatch_gaze.readthedocs.io/en/latest/api.html


=======
support
=======

All bugs, concerns and enhancement requests for this software can be submitted
here_.
All contributions, be it in the form of an issue or a pull-request,
are always welcome.


.. _here: https://github.com/adswa/multimatch_gaze/issues/new


================
acknowledgements
================

We thank Dr. Richard Dewhurst for kindly and swiftly providing the original
Matlab code for the MultiMatch toolbox via e-mail and being supportive of an
open source implementation.
