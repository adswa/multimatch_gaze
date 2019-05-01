import numpy as np
import pandas as pd
import os.path
import random
import collections
from bisect import bisect_right
from bisect import bisect_left
from .. import multimatch_gaze as mp

dtype = [('onset', '<f8'), ('duration', '<f8'),
         ('label', '<U10'), ('start_x', '<f8'),
         ('start_y', '<f8'), ('end_x', '<f8'),
         ('end_y', '<f8'), ('amp', '<f8'),
         ('peak_vel', '<f8'), ('med_vel', '<f8'), ('avg_vel', '<f8')]

def same_sample(run=1, subj=1):
    """duplicate dataset to force exactly similar scanpaths. Choose the run
    (integer between 1-8) and whether you want a lab (1) or mri (2) subject"""
    if subj == 1:
        sub = "sub-30"
    else:
        sub = "sub-10"
    path = os.path.join("multimatch_gaze/tests/testdata",
                        "{}_task-movie_run-{}_events.tsv".format(sub, run))
    loc = os.path.join("multimatch_gaze/tests/testdata",
                       "locations_run-{}_events.tsv".format(run))
    data = np.recfromcsv(path,
                         delimiter='\t',
                         dtype={'names': ('onset', 'duration', 'label',
                                          'start_x', 'start_y', 'end_x',
                                          'end_y', 'amp', 'peak_vel',
                                          'med_vel', 'avg_vel'),
                                'formats': ('f8', 'f8', 'U10', 'f8', 'f8',
                                            'f8', 'f8', 'f8', 'f8', 'f8',
                                            'f8')})
    data2 = data
    shots = pd.read_csv(loc, sep='\t')
    return data, data2, shots


def short_shots(run=3):
    """create a shortened shots location annotation to test longshots()"""
    loc = os.path.join("multimatch_gaze/tests/testdata",
                       "locations_run-{}_events.tsv".format(run))
    shots = pd.read_csv(loc, sep='\t')
    shortshots = shots[0:20]
    return shortshots


def mk_fix_vector(length=5):
    """creates a random length x 3 fixation vector in form of a record array"""
    fix = np.recarray((0,), dtype=[('start_x', '<f8'), ('start_y', '<f8'),
                                   ('duration', '<f8')])
    for i in range(0, length):
        fixation = np.array((np.random.uniform(1, 720),
                             np.random.uniform(1, 720),
                             np.random.uniform(0.01, 5)),
                            dtype=[('start_x', float),
                                   ('start_y', float),
                                   ('duration', float)])
        fix = np.append(fix, fixation)
    return fix


def mk_strucarray(length=5):
    """create a random scanpath in the data format generateScanpathStructureArray
    would output"""
    fixation_x = random.sample(range(700), length)
    fixation_y = random.sample(range(700), length)
    fixation_dur = random.sample(range(5), length)
    saccade_x = random.sample(range(700), length - 1)
    saccade_y = random.sample(range(700), length - 1)
    saccade_lenx = random.sample(range(700), length - 1)
    saccade_leny = random.sample(range(700), length - 1)
    saccade_rho = random.sample(range(700), length - 1)
    saccade_theta = random.sample(range(4), length - 1)
    eyedata = dict(
        fix=dict(
            x=fixation_x,
            y=fixation_y,
            dur=fixation_dur,
        ),
        sac=dict(
            x=saccade_x,
            y=saccade_y,
            lenx=saccade_lenx,
            leny=saccade_leny,
            theta=saccade_theta,
            rho=saccade_rho,
        )
    )
    eyedata2 = dict(
        fix=dict(
            x=fixation_x[::-1] * 2,
            y=fixation_y[::-1] * 2,
            dur=fixation_dur[::-1] * 2,
        ),
        sac=dict(
            x=saccade_x[::-1] * 2,
            y=saccade_y[::-1] * 2,
            lenx=saccade_lenx[::-1] * 2,
            leny=saccade_leny[::-1] * 2,
            theta=saccade_theta[::-1] * 2,
            rho=saccade_rho[::-1] * 2,
        )
    )
    return eyedata, eyedata2


def mk_angles():
    """creates vectors with predefined angular relations. angles1 and angles2
    contain the following properties: 1. same 0, 2. 60 diff, 3. 90 diff,
    4.120 diff,4. 180 diff (max. dissimilar). They are in sectors (0,1) and
    (0, -1).
    Angles3 and angles4 contain the same properties reversed and lie in sectors
    (-1, 0) and (-1, -1)"""
    angles1 = dict(sac=dict(theta=[0, 0.523, 0.785, 1.04, 1.57]))
    angles2 = dict(sac=dict(theta=[0, -0.523, -0.785, -1.04, -1.57]))
    angles3 = dict(sac=dict(theta=[1.57, 2.093, 2.356, 2.617, 3.14]))
    angles4 = dict(sac=dict(theta=[-1.57, -2.093, -2.356, -2.617, -3.14]))
    path = [0, 6, 12, 18, 24]
    M_assignment = np.arange(5 * 5).reshape(5, 5)
    return M_assignment, path, angles1, angles2, angles3, angles4


def mk_durs():
    """create some example duration for test_durationsim()"""
    durations1 = collections.OrderedDict()
    durations2 = collections.OrderedDict()
    durations1 = dict(fix=dict(dur=[0.001, 20.0, 7, -18, -2.0]))
    durations2 = dict(fix=dict(dur=[0.008, 18.0, 7, -11, 3.0]))
    path = [0, 6, 12, 18, 24]
    M_assignment = np.arange(5 * 5).reshape(5, 5)
    return M_assignment, path, durations1, durations2


def mk_supershort_shots():
    data = {'onset': np.arange(0, 20), 'duration': np.repeat(1, 20),
            'locale': np.repeat('somewhere', 20)}
    shots = pd.DataFrame(data)
    return shots


def mk_longershots():
    data = {'onset': np.arange(0, 20), 'duration': np.repeat(5, 20),
            'locale': np.repeat('somewhere', 20)}
    shots = pd.DataFrame(data)
    return shots


# some functions to work specifically with studyforrest eye tracking data

# Functions specifically for the data at hand

def takeclosestright(mylist, mynumber):
    """Return integer closest right to 'myNumber' in an ordered list.

    :param: mylist: int
    :param: mynumber: array

    :return: after: float, number within mylist closest to right of my number

    """

    pos = bisect_right(mylist, mynumber)
    if pos == 0:
        return mylist[0]
    if pos == len(mylist):
        return mylist[-1]
    after = mylist[pos]
    return after


def takeclosestleft(mylist, mynumber):
    """Return integer closest left to 'myNumber' in an ordered list.

    :param: mylist: int
    :param: mynumber: array

    :return: after: float, number within mylist closest to the left of mynumber
    """

    pos = bisect_left(mylist, mynumber)
    if pos == 0:
        return mylist[0]
    if pos == len(mylist):
        return mylist[-1]
    before = mylist[pos - 1]
    return before


def create_onsets(data, dur):
    """Create shot onsets from studyforrests location annotation.

    Create onset times of all shots of at least 'dur' seconds of length.

    :param: data: dataframe
        location annotation from studyforrest
    :param: dur: float
        time in seconds a shot should at least be long

    :return: onsets: array-like, list of shot onset times
    """
    onsets = []
    for index, row in data.iterrows():
        if row['duration'] >= dur:
            onsets.append(row['onset'])
    return onsets


def create_offsets(data, dur):
    """Create shot offsets from studyforrests location annotation.

    Create offset times of all shots of at least 'dur' seconds of length


    :param: data: dataframe, location annotation from studyforrest
    :param: dur: float, time in seconds a shot should at least be long

    :return: onsets: array-like, list of shot offset times
    """

    offsets = []
    for index, row in data.iterrows():
        if row['duration'] >= dur:
            # calculate end of shot by adding onset + duration, subtract an
            # epsilon to be really sure not to get into a cut
            offsets.append(row['onset'] + row['duration'] - 0.03)
    return offsets


def create_chunks(onsets, fixations, dur):
    """Chunk eyetracking data into scanpaths.

    Use onset data to obtain indices of full eyetracking data
    for chunking.

    :param: onsets:  array-like, onset times of movie shots
    :param: fixations: record array, nx4 fixation vector
        (onset, x, y, duration),
        output of preprocess() function
    :param: dur: float, desired duration of segment length

    :return: startidx, endix: array, start and end ids of eyemovement data
            to chunk into segments
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


def create_offsetchunks(offsets, fixations, dur):
    """Chunk eyetracking data into scanpaths.

    Use offset data to obtain indices of full eyetracking data
    for chunking.

    :param: offsets:  array-like, offset times of movie shots
    :param: fixations: record array, nx4 fixation vector
        (onset, x, y, duration), output of preprocess()
    :param: dur: float, desired duration of segment length

    :return: startidx, endix: array start and end ids of eyemovement data
            to chunk into segments
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
    """Chunk eyemovement data into scanpaths.

    :param: fixations: record array, nx4 fixation vector
        (onset, x, y, duration), output of preprocess()
    :param: startid, endid: array, start- and end-ids of the
        scanpaths, output from either create_chunks()
        or create_offsetchunks()

    :return: fixation_vector: array-like, a nx3 fixation vector
        (x, y, duration)
    """

    fixation_vector = []
    # slice fixation data according to indices, take columns
    # start_x, start_y and duration
    for idx in range(0, len(startid)):
        ind = fixations[startid[idx]:endid[idx]][["start_x",
                                                  "start_y",
                                                  "duration"]]
        fixation_vector.append(ind)
    return fixation_vector


def pursuits_to_fixations(remodnav_data):
    """Transform start and endpoints of pursuits to fixations.

    Uses the output of a record array created by the remodnav algorithm for
    eye-movement classification to transform pursuit data into fixations.
    The start and end point of a pursuit are relabeled as a fixation.
    This is useful for example if the underlying stimulus material is a
    moving image - visual intake of a moving object would then resemble
    a pursuit.

    :param: npdata: recordarray, remodnav output of eyemovement data

    :return: newdata: recordarray
    """
    # initialize empty rec array of the same shape
    newdata = np.recarray((0,), dtype=dtype)
    # reassemble rec array.
    # split pursuits to use end and start as fixations later
    from copy import deepcopy
    data = deepcopy(remodnav_data)
    for i, d in enumerate(data):
        if data[i]['label'] == 'PURS':
            # start and end point of pursuit get
            #  half the total duration
            d['duration'] = d['duration'] / 2
            d['label'] = 'FIXA'
            d2 = deepcopy(d)
            # end point of the pursuit is start
            # of new fixation
            d2['onset'] += d2['duration']
            d2['start_x'] = d2['end_x']
            d2['start_y'] = d2['end_y']
            newdata = np.append(newdata, np.array(d, dtype=dtype))
            newdata = np.append(newdata, np.array(d2, dtype=dtype))
        else:
            newdata = np.append(newdata, np.array(d, dtype=dtype))
    return newdata


def preprocess_remodnav(data, screensize):
    """Preprocess record array of eye-events.

    A record array from REMoDNaV data is preprocessed
    in the following way: Subset to only get fixation data,
    disregard out-of-frame gazes, subset to only keep x, y coordinates,
    duration.

    :param: data: recordarray, REMoDNaV output of eye events from movie
        data
    :param: screensize: list of float, screen measurements in px

    :return: fixations: array-like nx3 fixation vectors (onset, x, y,
        duration)

    """
    # only fixation labels
    filterevents = data[(data['label'] == 'FIXA')]
    # within x coordinates?
    filterxbounds = filterevents[np.logical_and(filterevents['start_x'] >= 0,
                                                filterevents['start_x'] <= screensize[0])]
    # within y coordinates?
    filterybounds = filterxbounds[np.logical_and(filterxbounds['start_y'] >= 0,
                                                 filterxbounds['end_y'] <= screensize[1])]
    # give me onset times, start_x, start_y and duration
    fixations = filterybounds[["onset", "start_x", "start_y",
                               "duration"]]
    return fixations

def read_remodnav(data):
    """ Helper to read input data produced by the REMoDNaV algorithm.
    Further information on the REMoDNaV algorithm can be found here:
    https://github.com/psychoinformatics-de/remodnav
    """
    d = np.recfromcsv(data,
        delimiter='\t',
        dtype=dtype
         )

    return d

def longshot(shots,
             group_shots,
             ldur=4.92):
    """Group movie shots without a cut together to obtain longer segments.

    Note: This way, fewer but longer scanpaths are obtained. Example: use
    median shotlength of 4.92s.

    :param: shots: dataframe, contains movie location annotation
    :param: group_shots: boolean, if True, grouping of movie shots is performed
    :param: dur: float, length in seconds for movie shot. An attempt is made to
        group short shots without a cut together to form longer shots of ldur
        length
    :return: aggregated, dataframe of aggregated movie shots
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


def docomparison_forrest(shots,
                         data1,
                         data2,
                         screensize=[1280, 720],
                         dur=4.92,
                         ldur=0,
                         offset=False,
                         TDur=0,
                         TDir=0,
                         TAmp=0,
                         grouping=False
                         ):
    """Compare two scanpaths on five similarity dimensions.

    :param: data1, data2: recarray, eyemovement information of forrest gump studyforrest dataset
    :param: screensize: list, screen dimensions in px.
    :param: ldur: float, duration in seconds. An attempt is made to group short shots
        together to form shots of ldur length
    :param: grouping: boolean, if True, simplification is performed based on thresholds TAmp,
        TDir, and TDur
    :param: TDir: float, Direction threshold, angle in degrees.
    :param: TDur: float, Duration threshold, duration in seconds.
    :param: TAmp: float, Amplitude threshold, length in px.

    :return: scanpathcomparisons: array
        array of 5 scanpath similarity measures
    :return: durations: array-like
        durations of extracted scanpaths. Vector (Shape), Direction
        (Angle), Length, Position, and Duration. 1 = absolute
        similarity, 0 = lowest similarity possible.
    :return: onsets: array-like
        onset times of the scanpaths
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
    fixations1 = preprocess_remodnav(newdata1, screensize)
    fixations2 = preprocess_remodnav(newdata2, screensize)
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
        exact_duration = fixations1[endid1[i]]['onset'] - \
                         fixations1[startid1[i]]['onset']
        # capture negative durations for invalid scanpaths
        if exact_duration > 0:
            exact_durations.append(exact_duration)
        else:
            exact_durations.append(np.nan)
        if i == len(startid1):
            print('Captured onsets and duration'
                  ' times of all scanpath pairs.')
    # loop over all fixation vectors/scanpaths and calculate similarity
    for i in range(0, len(onset)):
        # check if fixation vectors/scanpaths are long enough
        if (len(fixation_vectors1[i]) >= 3) &\
                (len(fixation_vectors2[i]) >= 3):
            subj1 = mp.gen_scanpath_structure(fixation_vectors1[i])
            subj2 = mp.gen_scanpath_structure(fixation_vectors2[i])
            if grouping:
                subj1 = mp.simplify_scanpath(subj1, TAmp, TDir, TDur)
                subj2 = mp.simplify_scanpath(subj2, TAmp, TDir, TDur)
            M = mp.cal_vectordifferences(subj1, subj2)
            scanpath_dim = np.shape(M)
            M_assignment = np.arange(scanpath_dim[0] *
                                     scanpath_dim[1]).reshape(scanpath_dim[0], scanpath_dim[1])
            weightedGraph = mp.createdirectedgraph(scanpath_dim, M, M_assignment)
            path, dist = mp.dijkstra(weightedGraph, 0, scanpath_dim[0] * scanpath_dim[1] - 1)
            unnormalised = mp.getunnormalised(subj1, subj2, path, M_assignment)
            normal = mp.normaliseresults(unnormalised, screensize)
            scanpathcomparisons.append(normal)
        # return nan as result if at least one scanpath it too short
        else:
            scanpathcomparisons.append(np.repeat(np.nan, 5))
            print('Scanpath {} had a length of {}, however, a minimal '
                  'length of 3 is required. Appending nan.'.format(i,
                                                                   min(len(fixation_vectors1), len(fixation_vectors2))))
    return scanpathcomparisons, onset_times, exact_durations
