import sys, os
sys.path.insert(0, os.path.abspath('./'))
import numpy as np
from . import utils as ut
from .. import MultiMatch as M
from .. import MultiMatch_pure as Mp


def test_same_real_data_forrest(run = 1, subj = 1):
    """
    Tests the studyforrest-specific functionality: are two identical scanpaths
    from the movie identical in all scanpath dimensions?
    :param run: specify the run (choose between 1 - 8)
    :param subj: specify the subject (example data of lab subject (1) or mri subject (2) available)
    """
    data1, data2, shots = ut.same_sample(run, subj)
    segments, onset, duration = M.doComparison(shots, data1, data2)
    segmentfinal = np.array(segments)
    assert np.all(segmentfinal.all(1))

# TODO: Full test of MultiMatch_pure
def test_same_real_data():
    """
    Tests the pure MultiMatch functionality: are two identical scanpaths supplied as fixation
    vectors identical in all scanpath dimensions?
    """
    import os
    testfile = os.path.abspath('MultiMatch/tests/testdata/segment_5_sub-19.tsv')
    data1 = np.recfromcsv(testfile,
                          delimiter='\t',
                          dtype={'names':('start_x', 'start_y', 'duration'),
                                 'formats':('f8', 'f8', 'f8')})
    grouping=False
    kwargs={"grouping": False}
    # TODO: ask Kyle why this fails
    results = Mp.doComparison(data1, data1, **kwargs)
    assert np.all(results.all(1))

def test_StructureGeneration(length = 5):
    """
    Tests whether the fixation vector is transformed into data of expected dimensions.
    :param length: Specify the length of the fixation vector to be transformed
    """
    fix_vector = ut.mk_fix_vector(length)
    results = Mp.generateStructureArrayScanpath(fix_vector)
    assert len(results[0]) == len(fix_vector)
    assert len(results[1]) == len(fix_vector)
    assert len(results[2]) == len(fix_vector)
    assert len(results[3]) == len(fix_vector)-1
    assert len(results[4]) == len(fix_vector)-1
    assert len(results[5]) == len(fix_vector)-1
    assert len(results[6]) == len(fix_vector)-1
    assert len(results[7]) == len(fix_vector)-1
    assert len(results[8]) == len(fix_vector)-1

def test_calVectorDifferences(length = 5):
    """
    Tests whether vector difference calculation results in matrix of expected dimensions.
    :param length: specify the length of the structured array to be used
    """
    data1, data2 = ut.mk_strucArray(length=5)
    Matrix = Mp.calVectordifferences(data1, data2)
    assert Matrix.shape[0] == length -1
    assert Matrix.shape[1] == 2*length -2

def test_anglesim():
    """
    Tests whether predefined angular relations (quarters of the unit circle) yield the expected similarity results.
    """
    M_assignment, path, angles1, angles2, angles3, angles4 = ut.mk_angles()
    a1a2 = Mp.calAngularDifference(angles1, angles2, path, M_assignment)
    a3a4 = Mp.calAngularDifference(angles3, angles4, path, M_assignment)
    assert a1a2 == [0, 1.0459999999999994, 1.5700000000000003, 2.08,
    3.1400000000000006]
    assert a3a4 == [3.1400000000000006, 2.0971853071795863, 1.5711853071795865,
    1.0491853071795862, 0.0031853071795859833]


def test_durationsim():
    """
    Tests whether similar but not identical durations result in duration similarities greater zero
    """
    M_assignment, path, duration1, duration2 = ut.mk_durs()
    res = Mp.calDurationDifference(duration1, duration2, path, M_assignment)
    assert all(np.asarray(res)>= 0)

#TODO: position similarity test, vector length test, vector shape test

def test_closestright():
    """
    Tests whether closestright function works as intended and does not
    unexpectedly return a number left to the reference.
    """
    mylist = np.arange(0, 25)
    for i in range(0, len(mylist)):
        res = M.takeClosestright(mylist, i)
        if i == 0:
            assert res == 1
        elif i == len(mylist)-1:
            assert res == 24
        else:
            assert res == i+1


def test_closestleft():
    """
    Tests whether closestright function works as intended and does not
    unexpectedly return a number right to the reference.
    """
    mylist = np.arange(0, 25)
    for i in range(0, len(mylist)):
        res = M.takeClosestleft(mylist, i)
        if i == 0:
            assert res == 0
        elif i == len(mylist)-1:
            assert res == 23
        else:
            assert res == i-1


def test_createChunks(run=1, subj=1):
    """
    Tests chunking of studyforrest data into scenes based on shot annotation.
    Are start and end ids of shots the same length?
    Does no endid preceed a start id?
    :param run: specify the run (choose between 1 - 8)
    :param subj: specify the subject (example data of lab subject (1) or mri subject (2) available)
    """
    data1, data2, shots = ut.same_sample(run, subj)
    onsets = M.createOnsets(shots, 3)
    fixations = M.preprocess(data1)
    startid, endid = M.createChunks(onsets, fixations, 3)
    assert len(startid) == len(endid)
    trues = []
    for i in range(0, len(startid)):
        trues.append(startid[i]<=endid[i])
    assert all(trues)


def test_longshot():
    """
    Test whether longshots combines scenes of longer length into the correct number of shots.
    :return:
    """
    shots = ut.short_shots()
    newshots = M.longshot(shots, dur = 4.92)
    assert len(newshots) == len(shots)-2


