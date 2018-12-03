import numpy as np
import pandas as pd
import os.path
import random


def same_sample(run=1, subj=1):
    '''duplicate dataset to force exactly similar scanpaths. Choose the run
    (integer between 1-8) and whether you want a lab (1) or mri (2) subject'''
    if subj == 1:
        sub = "sub-30"
    else:
        sub = "sub-10"
    path = os.path.join("MultiMatch/tests/testdata",
                        "{}_task-movie_run-{}_events.tsv".format(sub, run))
    loc = os.path.join("MultiMatch/tests/testdata",
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
    '''create a shortened shots location annotation to test longshots()'''
    loc = os.path.join("MultiMatch/tests/testdata",
                       "locations_run-{}_events.tsv".format(run))
    shots = pd.read_csv(loc, sep='\t')
    shortshots = shots[0:20]
    return shortshots


def mk_fix_vector(length=5):
    '''creates a random length x 3 fixation vector in form of a record array'''
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


def mk_strucArray(length=5):
    '''create a random scanpath in the data format generateScanpathStructureArray
    would output'''
    fixation_x = random.sample(range(700), length)
    fixation_y = random.sample(range(700), length)
    fixation_dur = random.sample(range(5), length)
    saccade_x = random.sample(range(700), length - 1)
    saccade_y = random.sample(range(700), length - 1)
    saccade_lenx = random.sample(range(700), length - 1)
    saccade_leny = random.sample(range(700), length - 1)
    saccade_rho = random.sample(range(700), length - 1)
    saccade_theta = random.sample(range(4), length - 1)
    eyedata = [fixation_x, fixation_y, fixation_dur, saccade_x, saccade_y,
               saccade_lenx, saccade_leny, saccade_theta, saccade_rho]
    eyedata2 = [fixation_x[::-1] * 2, fixation_y[::-1] * 2, fixation_dur[::-1] * 2,
                saccade_x[::-1] * 2, saccade_y[::-1] * 2, saccade_lenx[::-1] * 2, saccade_leny[::-1] * 2,
                saccade_theta[::-1] * 2, saccade_rho[::-1] * 2]
    return eyedata, eyedata2


def mk_angles():
    """creates vectors with predefined angular relations. angles1 and angles2
    contain the following properties: 1. same 0, 2. 60 diff, 3. 90 diff,
    4.120 diff,4. 180 diff (max. dissimilar). They are in sectors (0,1) and
    (0, -1).
    Angles3 and angles4 contain the same properties reversed and lie in sectors
    (-1, 0) and (-1, -1)"""
    angles1 = [[], [], [], [], [], [], [], [0, 0.523, 0.785, 1.04, 1.57], []]
    angles2 = [[], [], [], [], [], [], [], [0, -0.523, -0.785, -1.04, -1.57], []]
    angles3 = [[], [], [], [], [], [], [], [1.57, 2.093, 2.356, 2.617, 3.14], []]
    angles4 = [[], [], [], [], [], [], [], [-1.57, -2.093, -2.356, -2.617, -3.14], []]
    path = [0, 6, 12, 18, 24]
    M_assignment = np.arange(5 * 5).reshape(5, 5)
    return M_assignment, path, angles1, angles2, angles3, angles4


def mk_durs():
    '''create some example duration for test_durationsim()'''
    durations1 = [[], [], [0.001, 20.0, 7, -18, -2.0], [], [], [], [], [], []]
    durations2 = [[], [], [0.008, 18.0, 7, -11, 3.0], [], [], [], [], [], []]
    path = [0, 6, 12, 18, 24]
    M_assignment = np.arange(5 * 5).reshape(5, 5)
    return (M_assignment, path, durations1, durations2)


def mk_supershort_shots():
    data = {'onset': np.arange(0, 20), 'duration': np.repeat(1, 20), 'locale': np.repeat('somewhere', 20)}
    shots=pd.DataFrame(data)
    return shots


def mk_longershots():
    data = {'onset': np.arange(0, 20), 'duration': np.repeat(5, 20), 'locale': np.repeat('somewhere', 20)}
    shots=pd.DataFrame(data)
    return shots
