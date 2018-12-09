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
    orcid: ??
    affiliation: 1
affiliations:
  - name: Psychoinformatics Lab, Institute of Psychology,
    Otto-von-Guericke Universität Magdeburg
  - index: 1
date: 15.12.2018
bibliography: paper.bib
---

# Summary

The similarity of scanpaths, the trace of eye-movements
in space and time, offers insights into commonalities
and differences of viewing behaviour within and between
participants. In addition to quantification of position
and order of eye-movements, a comparison between them
adds a useful dimension to the traditional analysis of
eyetracking data.

``multimatch in Python`` is a Python based
reimplementation of the MultiMatch toolbox for scanpath
comparison, originally developed by [@Jarodzka] and
implemented by [@Dewhurst] in Matlab.
The multimatch method represents scanpaths as geometrical
vectors in a two-dimensional space: Any scanpath is build
up of a vector sequence in which the vectors represent
saccades, and the start and end position of saccadic
vectors represent fixations. Two such vector sequences
are compared on the five dimensions “vector shape”, “vector
length (amplitude)”, “vector position”, “vector direction”,
and “fixation duration” for a multidimensional similarity
evaluation.

The implementation aims at providing an open source
alternative to the existing MultiMatch toolbox. It further
contains optional functionality to compute scanpath
similarities in the eyetracking data of the studyforrest
dataset [@studyforrest].

The software is available as a github repository [@multimatch].

# Acknowledgements

We thank Dr. Richard Dewhurst for kindly and swiftly providing
the original matlab code for the MultiMatch toolbox via e-mail
and being supportive of an open source implementation.

# References

