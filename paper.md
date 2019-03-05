---
title: 'MultiMatch in Python'
tags:
    - eyetracking
    - scanpath
    - fixation
    - saccade
    - Python
authors:
  - name: Adina Wagner
    orcid: 0000-0003-2917-3450
    affiliation: 1
  - name: TODO: who else??
  - name: Michael Hanke
    orcid: 0000-0001-6398-6370
    affiliation: 1
affiliations:
  - name: Institute of Neuroscience and Medicine, Brain & Behaviour (INM-7),
    Research Centre Jülich
  - index: 1
date: 15.04.2019
bibliography: paper.bib
---

# Summary

The similarity of scanpaths, the trace of eye-movements
in space and time, offers insights into commonalities
and differences of viewing behavior within and between
participants. In addition to the quantification of position
and order of a series of eye-movements, a comparison
between them adds a useful dimension to the traditional
analysis of eyetracking data. For example, scanpath
comparisons are used to study analogy-making [@french],
visual exploration and imagery [@Johansson], habituation
in repetitive visual search [@burmester], or spatial
attention allocation in dynamic scenes [@mital]. The method
is applied within individuals as a measure of change [@burmester],
or across samples to study group differences [@french].
Therefore, in recent years, interest
in the study eye movement sequences has sparked the development
of novel methodologies and algorithms to perform scanpath
comparisons. However, many of the contemporary scanpath
comparison algorithms are implemented in closed-source,
non-free software such as Matlab.

``multimatch`` is a Python based
reimplementation of the MultiMatch toolbox for scanpath
comparison, originally developed by [@Jarodzka] and
implemented by [@Dewhurst] in Matlab.
This algorithm represents scanpaths as geometrical
vectors in a two-dimensional space: Any scanpath is build
up of a vector sequence in which the vectors represent
saccades, and the start and end position of saccadic
vectors represent fixations. Two such vector sequences
are - after optional simplification based on angular relations
and amplitudes of saccades -
compared on the five dimensions “vector shape”, “vector
length (amplitude)”, “vector position”, “vector direction”,
and “fixation duration” for a multidimensional similarity
evaluation.

This reimplementation in Python aims at providing an
accessible, well-documented, and tested open
source alternative to the existing MultiMatch toolbox. The algorithm
is an established tool for scanpath comparison [@anderson],
and improved availability will likely foster its application
in a broader research community. multimatch
is available in the form of the github repository [@multimatch]
and as the python module multimatch via 'pip install multimatch'.
The module contains the same functionality as the original
Matlab toolbox, that is, scanpath comparison with optional
simplification according to user-defined thresholds, and it
provides this functionality within command-line calls or
interactive python sessions.

# Acknowledgements

We thank Dr. Richard Dewhurst for kindly and swiftly providing
the original matlab code for the MultiMatch toolbox via e-mail
and being supportive of an open source implementation.

# References

