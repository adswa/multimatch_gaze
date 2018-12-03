# emacs: -*- mode: python; py-indent-offset: 4; tab-width: 4; indent-tabs-mode: nil -*-
# ex: set sts=4 ts=4 sw=4 noet:
# ## ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the multimatch package for the
#   copyright and license terms.
#
# ## ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

from __future__ import absolute_import
from .multimatch import docomparison
from .multimatch_forrest import doComparisonForrest

__version__ = '0.1'


#
# import logging
# lgr = logging.getLogger('multimatch')
#
# def multimatch_pure():#args=sys.argv):
#     import argparse
#     import inspect
#
#     parser = argparse.ArgumentParser()
#     parser.add_argument(
#         'input1', metavar='<datafile>',
#         help= """eyemovement data of the first scanpath. The first
#         two columns in this file must contain x and y coordinates, the third
#         column must contain durations. The file should be tab seperated. This
#         file is read with NumPys recfromcsv.""")
#     parser.add_argument(
#         'input2', metavar='<datafile>',
#         help= """eyemovement data of the second scanpath. The first
#         two columns in this file must contain x and y coordinates, the third
#         column must contain durations. The file should be tab seperated. This
#         file is read with NumPys recfromcsv.""")
#     parser.add_argument(
#         'direction_threshold', type=float, default=0.0, metavar='<TDir>',
#         help="""Direction threshold for direction based grouping in degrees.
#         If 0, no simplification will be performed."""
#     )
#     parser.add_argument(
#         'amplitude_threshold', type=float, default=0.0, metavar='<TAmp>',
#         help="""Amplitude threshold for length based grouping in pixel.
#         If 0, no simplification will be performed."""
#     )
#     parser.add_argument(
#         'duration_threshold', type=float, default=0.0, metavar='<TDur>',
#         help="""Duration threshold for scanpath simplification."""
#     )
#     parser.add_argument(
#         'screensize', default=[1280, 720], metavar='<sz>',
#         help="""Resolution of screen in px, default is [1280, 720]."""
#     )
#
#     args = parser.parse_args()
#
#     #read in data
#     data1 = np.recfromcsv(
#         args.input1,
#         delimiter='\t',
#         dtype={'names':('start_x', 'start_y', 'duration'),
#                'formats':('f8', 'f8', 'f8')})
#     data2 = np.recfromcsv(
#         args.input2,
#         delimiter='\t',
#         dtype={'names':('start_x', 'start_y', 'duration'),
#                'formats':('f8', 'f8', 'f8')})
#
#     TDir = args.direction_threshold
#     TAmp = args.amplitude_threshold
#     TDur = args.duration_threshold
#     sz = args.screensize
#
#     if (TDir != 0) and (TAmp != 0):
#         grouping = True
#         print('Scanpath comparison is done with grouping saccades shorter than {}px and with an angle smaller than {} degrees'
#               ' if consecutive fixation are shorter than {} seconds.'.format(TAmp, TDir, TDur))
#     else:
#         grouping = False
#         print('Scanpath comparison is done without any grouping')
#
# if __name__ == '__main__':
#     main()