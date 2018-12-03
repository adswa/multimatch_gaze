#!/usr/bin/python3


import numpy as np
import pandas as pd
from bisect import bisect_right
from bisect import bisect_left
import sys
import os

sys.path.insert(0, os.path.abspath('./'))
import MultiMatch.multimatch as Mp


# Functions specifically for the data at hand

def takeclosestright(mylist, mynumber):
    """Return integer closest right to 'myNumber' in an ordered list

    Parameters
    ------------
    mylist: int
    mynumber: array

    Returns
    ------------
    after: float
        number within mylist closest to the right of mynumber

    References
    -----------

    I stole this from
    https://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-
    closest-to-a-given-value

    """

    pos = bisect_right(mylist, mynumber)
    if pos == 0:
        return mylist[0]
    if pos == len(mylist):
        return mylist[-1]
    after = mylist[pos]
    return after


def takeclosestleft(myList, myNumber):
    """Return integer closest left to 'myNumber' in an ordered list

    Parameters
    ------------
    mylist: int
    mynumber: array

    Returns
    ------------
    after: float
        number within mylist closest to the left of mynumber


    """

    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    return before


def create_onsets(data, dur):
    """Create shot onsets from studyforrests location annotation

    Create onset times of all shots of at least 'dur' seconds of length

    Parameters
    -----------

    data: dataframe
        location annotation from studyforrest
    dur: float
        time in seconds a shot should at least be long

    Returns
    ----------
    onsets: array-like
        list of shot onset times
    """
    onsets = []
    for index, row in data.iterrows():
        if row['duration'] >= dur:
            onsets.append(row['onset'])
    return onsets


def create_offsets(data, dur):
    """Create shot offsets from studyforrests location annotation

    Create offset times of all shots of at least 'dur' seconds of length

    Parameters
    -----------

    data: dataframe
        location annotation from studyforrest
    dur: float
        time in seconds a shot should at least be long

    Returns
    ----------
    onsets: array-like
        list of shot offset times
    """

    offsets = []
    for index, row in data.iterrows():
        if row['duration'] >= dur:
            # calculate end of shot by adding onset + duration, subtract an
            # epsilon to be really sure not to get into a cut
            offsets.append(row['onset'] + row['duration'] - 0.03)
    return offsets


def create_chunks(onsets, fixations, dur):
    """Chunk eyetracking data into scanpaths

    Use onset data to obtain indices of full eyetracking data
    for chunking.

    Parameters
    -----------

    onsets:  array-like
        onset times of movie shots
    fixations: record array
        nx4 fixation vector (onset, x, y, duration),
        output of preprocess() function
    dur: float
        desired duration of segment length

    Returns
    -----------
    startidx, endix: array
        start and end ids of eyemovement data to chunk into segments
    """
    # initialize empty lists
    startidx, endidx = [], []
    for shotonset in onsets:
        start = takeclosestright(fixations['onset'], shotonset)
        startidx.append(np.where(fixations['onset'] == start)[0].tolist())
        end = takeclosestright(fixations['onset'], shotonset + dur)
        endidx.append(np.where(fixations['onset'] == end)[0].tolist())
    # flatten the nested lists
    startidx = [element for sublist in startidx for element in sublist]
    endidx = [element for sublist in endidx for element in sublist]
    return startidx, endidx


# createOffsetChunks changes to **last** seconds of shot!!

def create_offsetchunks(offsets, fixations, dur):
    """Chunk eyetracking data into scanpaths

    Use offset data to obtain indices of full eyetracking data
    for chunking.

    Parameters
    -----------

    offsets:  array-like
        offset times of movie shots
    fixations: record array
        nx4 fixation vector (onset, x, y, duration), output of preprocess()
    dur: float
        desired duration of segment length

    Returns
    -----------
    startidx, endix: array
        start and end ids of eyemovement data to chunk into segments
    """
    startidx, endidx = [], []
    for shotoffset in offsets:
        start = takeclosestright(fixations['onset'], shotoffset - dur)
        startidx.append(np.where(fixations['onset'] == start)[0].tolist())
        end = takeclosestleft(fixations['onset'], shotoffset)
        endidx.append(np.where(fixations['onset'] == end)[0].tolist())
    # flatten the nested lists
    startidx = [element for sublist in startidx for element in sublist]
    endidx = [element for sublist in endidx for element in sublist]
    return startidx, endidx


def fixations_chunks(fixations, startid, endid):
    """Chunk eyemovement data into scanpaths

    Parameters
    -----------
    fixations: record array
        nx4 fixation vector (onset, x, y, duration), output of preprocess()
    startid, endid: array
        start- and end-ids of the scanpaths, output from either create_chunks()
        or create_offsetchunks()

    Returns
    -----------
    fixation_vector: array-like
        a nx3 fixation vector (x, y, duration)
    """

    fixation_vector = []
    # slice fixation data according to indices, take columns start_x, start_y and
    # duration
    for idx in range(0, len(startid)):
        ind = fixations[startid[idx]:endid[idx]][["start_x",
                                                  "start_y",
                                                  "duration"]]
        fixation_vector.append(ind)
    return fixation_vector


def pursuits_to_fixations(npdata):
    """Transform start and endpoints of pursuits to fixations

    Uses the output of a record array created by the remodnav algorithm for eye-
    movement classification to transform pursuit data into fixations.

    Parameters
    ------------
    npdata: recordarray
        remodnav output of eyemovement data

    Returns
    -----------
    newdata: recordarray
    """
    # initialize empty rec array of the same shape
    newdata = np.recarray((0,), dtype=[('onset', '<f8'),
                                       ('duration', '<f8'),
                                       ('label', '<U10'),
                                       ('start_x', '<f8'),
                                       ('start_y', '<f8'),
                                       ('end_x', '<f8'),
                                       ('end_y', '<f8'),
                                       ('amp', '<f8'),
                                       ('peak_vel', '<f8'),
                                       ('med_vel', '<f8'),
                                       ('avg_vel', '<f8')])
    # reassemble rec array.
    # split pursuits to use end and start as fixations later
    for i in range(0, len(npdata)):
        if npdata[i]['label'] == 'PURS':
            row_1 = npdata[i]
            row_1['duration'] = npdata[i]['duration'] / 2
            row_2 = row_1.copy()
            row_2['onset'] += row_2['duration']
            row_2['start_x'] = row_2['end_x']
            row_2['start_y'] = row_2['end_y']
            newdata = np.append(newdata, row_1)
            newdata = np.append(newdata, row_2)
        else:
            newdata = np.append(newdata, npdata[i])
    return newdata


def preprocess(data, sz=[1280, 720]):
    """Preprocess record array of eye-events

    A record array of the studyforrest eyemovement data is preprocessed
    in the following way: Subset to only get fixation and pursuit data,
    disregard out-of-frame gazes, subset to only keep x, y coordinates,
    duration, and onset.

    Parameters
    ------------
    data: recordarray
        remodnav output of eye events from movie data
    sz: list of float
        screen measurements in px

    Returns
    -----------
    fixations: array-like
        nx4 fixation vectors (onset, x, y, duration)

    """

    # only fixations and pursuits
    Filterevents = data[np.logical_or(data['label'] == 'FIXA',
                                      data['label'] == 'PURS')]
    # within x coordinates?
    Filterxbounds = Filterevents[np.logical_and(Filterevents['start_x'] >= 0,
                                                Filterevents['start_x'] <= sz[0])]
    # within y coordinates?
    Filterybounds = Filterxbounds[np.logical_and(Filterxbounds['start_y'] >= 0,
                                                 Filterxbounds['end_y'] <= sz[1])]
    # give me onset times, start_x, start_y and duration
    fixations = Filterybounds[["onset", "start_x", "start_y",
                               "duration"]]
    return fixations


def longshot(shots,
             group_shots,
             ldur=4.92):
    """Group movie shots without a cut together to obtain longer segments

    Note: This way, fewer but longer scanpaths are obtained. Example: use
    median shotlength of 4.92s.

    Parameters
    ------------
    shots: dataframe
        contains movie location annotation
    group_shots: boolean
        if True, grouping of movie shots is performed
    dur: float
        length in seconds for movie shot. An attempt is made to group short
        shots without a cut together to form longer shots of ldur length
    :return:
    """
    # turn pandas dataframe shots into record array
    structshots = shots.to_records()
    if group_shots:
        i = 0
        while i < len(structshots):
            # break before running into index error
            if structshots[i] == structshots[-1]:
                break
            else:
                if (structshots[i]['duration'] < ldur) & \
                        (structshots[i + 1]['duration'] < ldur) & \
                        (structshots[i]['locale'] == structshots[i + 1]['locale']):
                    # add durations together and delete second row
                    structshots[i]['duration'] += structshots[i + 1]['duration']
                    structshots = np.delete(structshots, i + 1, 0)
                else:
                    i += 1
    aggregated = pd.DataFrame({'onset': structshots['onset'].tolist(),
                               'duration': structshots['duration'].tolist()},
                              columns=['onset', 'duration'])
    return aggregated


def doComparisonForrest(shots,
                        data1,
                        data2,
                        sz=[1280, 720],
                        dur=5,
                        ldur=0.0,
                        offset=False,
                        TDur=0,
                        TDir=0,
                        TAmp=0,
                        grouping=False
                        ):
    """Compare two scanpaths on five similarity dimensions

     Parameters
     ------------

    data1, data2: recarray
        eyemovement information of forrest gump studyforrest dataset
    sz: list
        screen dimensions in px.
    ldur: float
        duration in seconds. An attempt is made to group short shots
        together to form shots of ldur length
    grouping: boolean
        if True, simplification is performed based on thresholds TAmp,
        TDir, and TDur
    TDir: float
        Direction threshold, angle in degrees.
    TDur: float
        Duration threshold, duration in seconds.
    TAmp: float
        Amplitude threshold, length in px.

    Returns
    ------------
    scanpathcomparisons: array
        array of 5 scanpath similarity measures
    durations: array-like
        durations of extracted scanpaths. Vector (Shape), Direction
        (Angle), Length, Position, and Duration. 1 = absolute
        similarity, 0 = lowest similarity possible.
    onsets: array-like
        onset times of the scanpaths

    Examples
    ------------
    """
    # determine whether short shots should be grouped together
    if ldur != 0:
        group_shots = True
    else:
        group_shots = False

    scanpathcomparisons = []
    # transform pursuits into fixations
    newdata1 = pursuits_to_fixations(data1)
    newdata2 = pursuits_to_fixations(data2)
    print('Loaded data.')
    # preprocess input files
    fixations1 = preprocess(newdata1, sz)
    fixations2 = preprocess(newdata2, sz)
    shots = longshot(shots, group_shots, ldur)
    # get shots and scanpath on- and offsets
    if offset:
        onset = create_offsets(shots, dur)
        startid1, endid1 = create_offsetchunks(onset, fixations1, dur)
        startid2, endid2 = create_offsetchunks(onset, fixations2, dur)
    else:
        onset = create_onsets(shots, dur)
        startid1, endid1 = create_chunks(onset, fixations1, dur)
        startid2, endid2 = create_chunks(onset, fixations2, dur)
    fixation_vectors1 = fixations_chunks(fixations1, startid1, endid1)
    fixation_vectors2 = fixations_chunks(fixations2, startid2, endid2)
    print('Split fixation data into {} scanpaths.'.format(len(startid1)))
    # save onset and duration times, if valid ones can be calculated
    onset_times = []
    exact_durations = []
    for i in range(0, len(startid1)):
        onset_time = fixations1[startid1[i]]['onset']
        onset_times.append(onset_time)
        exact_duration = fixations1[endid1[i]]['onset'] - fixations1[startid1[i]]['onset']
        # capture negative durations for invalid scanpaths
        if exact_duration > 0:
            exact_durations.append(exact_duration)
        else:
            exact_durations.append(np.nan)
        if i == len(startid1):
            print('Captured onsets and duration times of all scanpath pairs.')
    # loop over all fixation vectors/scanpaths and calculate similarity
    for i in range(0, len(onset)):
        # check if fixation vectors/scanpaths are long enough
        if (len(fixation_vectors1[i]) >= 3) & (len(fixation_vectors2[i]) >= 3):
            print('Computing similarity for comparison {}.'.format(i+1))
            subj1 = Mp.gen_scanpath_structure(fixation_vectors1[i])
            subj2 = Mp.gen_scanpath_structure(fixation_vectors2[i])
            if grouping:
                subj1 = Mp.simplify_scanpath(subj1, TAmp, TDir, TDur)
                subj2 = Mp.simplify_scanpath(subj2, TAmp, TDir, TDur)
                print('Simplification of pair {} completed.'.format(i))
            M = Mp.cal_vectordifferences(subj1, subj2)
            szM = np.shape(M)
            M_assignment = np.arange(szM[0] * szM[1]).reshape(szM[0], szM[1])
            weightedGraph = Mp.createdirectedgraph(szM, M, M_assignment)
            path, dist = Mp.dijkstra(weightedGraph, 0, szM[0] * szM[1] - 1)
            unnormalised = Mp.getunnormalised(subj1, subj2, path, M_assignment)
            normal = Mp.normaliseresults(unnormalised, sz)
            scanpathcomparisons.append(normal)
            print('Done.')
        # return nan as result if at least one scanpath it too short
        else:
            scanpathcomparisons.append(np.repeat(np.nan, 5))
            print('Scanpath {} had a length of {}, however, a minimal '
                  'length of 3 is required. Appending nan.'.format(i, min(len(fixation_vectors1), len(fixation_vectors2))))
    return scanpathcomparisons, onset_times, exact_durations

def main():
    print('Attempting to simplify scanpaths')
    segment, onset, duration = doComparisonForrest(shots,
                                                   data1,
                                                   data2,
                                                   sz,
                                                   dur,
                                                   ldur,
                                                   offset,
                                                   TDur,
                                                   TDir,
                                                   TAmp,
                                                   grouping)
    segmentfinal = np.array(segment)
    results = np.column_stack((onset, duration, segmentfinal))
    # save
    print('Saving results at {}.'.format(args.output))
    if not os.path.isdir(os.path.dirname(args.output)):
        os.makedirs(os.path.dirname(args.output))
    np.savetxt(args.output,
               results,
               fmt='%f\t%f\t%f\t%f\t%f\t%f\t%f',
               delimiter='\t',
               header="onset\tduration\tvector_sim\tdirection_sim\tlength_sim\tposition_sim\tduration_sim",
               comments=''
               )



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    # define arguments
    parser.add_argument('-i', '--input1', nargs='+', help='Input1: eyemovement data of the first subject',
                        metavar='PATH', required=True)
    parser.add_argument('-j', '--input2', nargs='+', help='Input2: eyemovement data of the second subject',
                        metavar='PATH', required=True)
    parser.add_argument('-k', '--input3', help='Input3: location annotation of the movie segment', metavar='PATH',
                        required=True)
    parser.add_argument('-o', '--output', help='Output: Specify path where output should be saved', metavar='PATH',
                        required=True)
    parser.add_argument('-d', '--duration',
                        help='duration: Specify the approx. time of a scanpath in seconds, i.e. 3.0. Note: Scanpaths are extracted within a shot, not across shots! Long durations will lead to only few scanpaths',
                        type=float, default=4.92)
    parser.add_argument('-ld', '--lduration',
                        help='duration: group short shots in the same locale (i.e. no change of scenes between them) together for longer scanpaths.',
                        type=float, default=0.0)
    parser.add_argument('-di', '--direction_threshold',
                        help='direction_threshold: for direction based grouping. If 0: no grouping will be performed',
                        type=float, default=0.0)
    parser.add_argument('-am', '--amplitude_threshold',
                        help='amplitude_threshold: for amplitude based grouping. If 0: no grouping will be performed',
                        type=float, default=0.0)
    parser.add_argument('-du', '--duration_threshold', help='duration_threshold: for direction based grouping.',
                        type=float, default=0.0)
    parser.add_argument('-sz', '--screensize', help='screensize: Resolution of screen in px, default is [1280, 720]',
                        default=[1280, 720])
    parser.add_argument('-pos', '--position_offset',
                        help='position_offset: if True, scanpaths of dur length stop at shotoffset (instead of starting at shotonset',
                        default=False)

    args = parser.parse_args()

    # Data read in

    shots = pd.read_csv(args.input3, sep='\t')
    data1 = np.recfromcsv(args.input1[0],
                          delimiter='\t',
                          dtype={'names': ('onset', 'duration', 'label', 'start_x', 'start_y',
                                           'end_x', 'end_y', 'amp', 'peak_vel', 'med_vel', 'avg_vel'),
                                 'formats': ('f8', 'f8', 'U10', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',
                                             'f8', 'f8')})

    data2 = np.recfromcsv(args.input2[0],
                          delimiter='\t',
                          dtype={'names': ('onset', 'duration', 'label', 'start_x', 'start_y',
                                           'end_x', 'end_y', 'amp', 'peak_vel', 'med_vel', 'avg_vel'),
                                 'formats': ('f8', 'f8', 'U10', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',
                                             'f8', 'f8')})
    dur = args.duration
    TDir = args.direction_threshold
    TAmp = args.amplitude_threshold
    TDur = args.duration_threshold
    sz = args.screensize
    ldur = args.lduration
    offset = args.position_offset
    # derive simple boolean variable to
    if (TDir != 0) and (TAmp != 0):
        grouping = True
        print(
            'Scanpath comparison is done with grouping saccades shorter than {}px and with an angle smaller than {} degrees'
            ' if consecutive fixation are shorter than {} seconds.'.format(TAmp, TDir, TDur))
    else:
        grouping = False
        print('Scanpath comparison is done without any grouping')

    # Execution
    main()

#    segment, onset, duration = doComparisonForrest(shots,
#                                                   data1,
#                                                   data2,
#                                                   sz,
#                                                   dur,
#                                                   ldur,
#                                                   offset,
#                                                   TDur,
#                                                   TDir,
#                                                   TAmp,
#                                                   grouping)
#    segmentfinal = np.array(segment)
#    results = np.column_stack((onset, duration, segmentfinal))
#    # save
#    if not os.path.isdir(os.path.dirname(args.output)):
#        os.makedirs(os.path.dirname(args.output))
#    np.savetxt(args.output,
#               results,
#               fmt='%f\t%f\t%f\t%f\t%f\t%f\t%f',
#               delimiter='\t',
#               header="onset\tduration\tvector_sim\tdirection_sim\tlength_sim\tposition_sim\tduration_sim",
#               comments=''
#               )
