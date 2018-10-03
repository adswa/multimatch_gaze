import sys, os
sys.path.insert(0, os.path.abspath('..'))
import pytest
import numpy as np
from . import utils as ut
from .. import MultiMatch as M

def test_same_real_data(run = 1, subj = 1):
    data1, data2, shots = ut.same_sample(run, subj)
    segments, onset, duration = M.doComparison(shots, data1, data2)
    segmentfinal = np.array(segments)
    assert np.all(segmentfinal.all(1))

def test_StructureGeneration(length = 5):
    fix_vector = ut.mk_fix_vector(length)
    results = M.generateStructureArrayScanpath(fix_vector)
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
    data1, data2 = ut.mk_strucArray(length=5)
    Matrix = M.calVectordifferences(data1, data2)
    assert Matrix.shape[0] == length -1
    assert Matrix.shape[1] == 2*length -2


def test_anglesim():
    M_assignment, path, angles1, angles2, angles3, angles4 = ut.mk_angles()
    a1a2 = M.calAngularDifference(angles1, angles2, path, M_assignment)
    a3a4 = M.calAngularDifference(angles3, angles4, path, M_assignment)
    assert a1a2 == [0, 1.0459999999999994, 1.5700000000000003, 2.08,
    3.1400000000000006]
    assert a3a4 == [3.1400000000000006, 2.0971853071795863, 1.5711853071795865,
    1.0491853071795862, 0.0031853071795859833]

def test_closestright():
    mylist = np.arange(0, 25)
    for i in range(0, len(mylist)):
        res = M.takeClosestright(mylist, i)
        if i == 0:
            assert res == 1
        elif i == len(mylist)-1:
            assert res == 24
        else:
            assert res == i+1


def test_createChunks(run=1, subj=1):
    data1, data2, shots = ut.same_sample(run, subj)
    onsets = M.createOnsets(shots, 3)
    fixations = M.preprocess(data1)
    startid, endid = M.createChunks(onsets, fixations, 3)
    assert len(startid) == len(endid)
    trues = []
    for i in range(0, len(startid)):
        trues.append(startid[i]<=endid[i])
    assert all(trues)


def test_durationsim():
    M_assignment, path, duration1, duration2 = ut.mk_durs()
    res = M.calDurationDifference(duration1, duration2, path, M_assignment)
    assert all(np.asarray(res)>= 0)


def test_longshot():
    shots = ut.short_shots()
    newshots = M.longshot(shots, dur = 4.92)
    assert len(newshots) == len(shots)-2


