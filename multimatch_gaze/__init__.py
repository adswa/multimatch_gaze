# emacs: -*- mode: python; py-indent-offset: 4; tab-width: 4; indent-tabs-mode: nil -*-
# ex: set sts=4 ts=4 sw=4 noet:
# ## ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See LICENSE file distributed along with the multimatch-gaze package for the
#   license terms.
#
# ## ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

# this import is necessary for the tests to run without having multimatch-gaze installed
# locally - else the modules can't be imported
from __future__ import absolute_import

__version__ = '0.1.1'

from .multimatch_gaze import docomparison
from .multimatch_gaze import remodnav_reader
