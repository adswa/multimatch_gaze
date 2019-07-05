---
title: 'multimatch_gaze: The MultiMatch algorithm for gaze path comparison in Python'
tags:
  - eyetracking
  - scan path
  - fixation
  - saccade
  - Python
authors:
 - name: Adina S. Wagner
   orcid: 0000-0003-2917-3450
   affiliation: 1
 - name: Yaroslav O. Halchenko
   orcid: 0000-0003-3456-2493
   affiliation: "3"
 - name: Michael Hanke
   orcid: 0000-0001-6398-6370
   affiliation: "1, 2"
affiliations:
 - name: Institute of Neuroscience and Medicine, Brain & Behaviour (INM-7), Research Centre Jülich, Jülich, Germany
   index: 1
 - name: Institute of Systems Neuroscience, Medical Faculty, Heinrich Heine University Düsseldorf, Düsseldorf, Germany
   index: 2
 - name: Department of Psychological and Brain Sciences, Dartmouth College, Dartmouth, NH, United States
   index: 3
date: 17 May 2019
bibliography: paper.bib
---

# Summary

The similarity of scan paths, the trace of eye-movements
in space and time, offers insights into commonalities
and differences of viewing behavior within and between
observers. In addition to the quantification of position
and order of a series of eye-movements, a comparison
between them adds an insightful dimension to the traditional
analysis of eyetracking data. For example, scan path
comparisons are used to study analogy-making [@french],
visual exploration and imagery [@Johansson], habituation
in repetitive visual search [@burmester], or spatial
attention allocation in dynamic scenes [@mital]. The method
is applied within individuals as a measure of change [@burmester],
or across samples to study group differences [@french].
Therefore, in recent years, interest
in the study of eye movement sequences has sparked the development
of novel methodologies and algorithms to perform scan path
comparisons. However, many of the contemporary scan path
comparison algorithms are implemented in closed-source,
non-free software such as Matlab.

``multimatch-gaze`` is a Python based
reimplementation of the MultiMatch toolbox for scan path
comparison, originally developed by [@Jarodzka] and
implemented by [@Dewhurst] in Matlab.
This algorithm represents scan paths as geometrical
vectors in a two-dimensional space: Any scan path is built
up of a coordinate vector sequence in which the start and end position
of vectors represent fixations, and the vectors represent
saccades. Two such vector sequences
are, after optional simplification based on angular relations
and amplitudes of saccades,
compared on the five dimensions “vector shape”, “vector
length (amplitude)”, “vector position”, “vector direction”,
and “fixation duration” for a multidimensional similarity
evaluation.

This reimplementation in Python aims at providing an
accessible, documented, and tested open
source alternative to the existing MultiMatch toolbox. The algorithm
is an established tool for scan path comparison [@anderson],
and improved availability aids adoption
in a broader research community. multimatch-gaze
is available from its Github repository
and as the Python package ``multimatch-gaze`` via ``pip install multimatch-gaze``.
The module contains the same functionality as the original
Matlab toolbox, that is, scan path comparison with optional
simplification according to user-defined thresholds, and it
provides this functionality via a command line interface or
a Python API.

Data for scan path comparison can be supplied as nx3
fixation vectors with columns corresponding to x-coordinates,
y-coordinates, and duration of the fixation in seconds (as for
the original Matlab toolbox).
Alternatively, multimatch-gaze can natively read in event detection
output produced by REMoDNaV [@remodnav], a velocity-based eye movement
classification algorithm written in Python.
For REMoDNaV-based input,
users can additionally specify whether smooth pursuit events
in the data should be kept in the scan path or
discarded.

# Acknowledgements

We thank Dr. Richard Dewhurst for kindly and swiftly providing
the original Matlab code for the MultiMatch toolbox via e-mail
and being supportive of an open source implementation.

# References

