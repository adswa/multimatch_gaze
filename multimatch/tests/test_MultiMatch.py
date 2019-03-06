import sys
import os
import numpy as np
from . import utils as ut
from .. import multimatch as mp


def test_same_real_data_forrest(run=1,
                                subj=1):
    """Test identical forrest scanpaths

    Tests whether two identical scanpaths
    from the movie identical in all scanpath dimensions(no grouping).

    Parameters
    -----------
    run: int
        specify the run (choose between 1 - 8)
    subj: int
        specify the subject (example data of lab subject (1)
        or mri subject (2) available)
    """
    data1, data2, shots = ut.same_sample(run, subj)
    segments, onset, duration = ut.docomparison_forrest(shots,
                                                        data1,
                                                        data2,
                                                        sz=[1280, 720],
                                                        dur=3,
                                                        ldur=0,
                                                        offset=False,
                                                        TDur=0.0,
                                                        TDir=0.0,
                                                        TAmp=0.0,
                                                        grouping=False
                                                        )
    segmentfinal = np.array(segments)
    assert np.all(segmentfinal.all(1))


def test_same_real_data():
    """Test identical scanpaths

    smoketest reading in of fixation vectors, tests whether two identical
    scanpaths supplied as fixation vectors identical in all scanpath
    dimensions?
    """
    import os
    testfile = os.path.abspath('multimatch/tests/testdata/segment_5_sub-19.tsv')
    data1 = np.recfromcsv(testfile,
                          delimiter='\t',
                          dtype={'names': ('start_x', 'start_y', 'duration'),
                                 'formats': ('f8', 'f8', 'f8')})
    results = mp.docomparison(data1,
                              data1,
                              sz=[720, 1280],
                              grouping=False,
                              TDir=0,
                              TDur=0,
                              TAmp=0)
    resultsfinal = np.array(results)
    assert np.all(resultsfinal.all(1))


def test_simplification(run=1,
                        subj=1):
    """
    Smoketest to see whether simplification blows up and whether identical
    scanpaths are still identical after grouping, using real data.
    """
    testfile = os.path.abspath('multimatch/tests/testdata/segment_5_sub-19.tsv')
    data1 = np.recfromcsv(testfile,
                          delimiter='\t',
                          dtype={'names': ('start_x', 'start_y', 'duration'),
                                 'formats': ('f8', 'f8', 'f8')})
    results = mp.docomparison(data1,
                              data1,
                              sz=[1280, 720],
                              grouping=True,
                              TDir=30.0,
                              TDur=0.05,
                              TAmp=100.0)
    Mdata1, Mdata2, shots = ut.same_sample(run, subj)
    segments, onset, duration = ut.docomparison_forrest(shots,
                                                        Mdata1,
                                                        Mdata2,
                                                        dur=3,
                                                        ldur=0.0,
                                                        offset=False,
                                                        sz=[1280, 720],
                                                        grouping=True,
                                                        TDir=30.0,
                                                        TDur=0.05,
                                                        TAmp=100.0)
    resultsfinal = np.array(results)
    resultsfinal_M = np.array(segments)
    assert np.all(resultsfinal.all(1))
    assert np.all(resultsfinal_M.all(1))


def test_structure_generation(length=5):
    """
    Tests whether the fixation vector is transformed into data of
    expected dimensions.

    Parameters
    -----------
    length: int
        length of the fixation vector to be transformed
    """
    fix_vector = ut.mk_fix_vector(length)
    results = mp.gen_scanpath_structure(fix_vector)
    assert len(results['fixation_x']) == len(fix_vector)
    assert len(results['fixation_y']) == len(fix_vector)
    assert len(results['fixation_dur']) == len(fix_vector)
    assert len(results['saccade_x']) == len(fix_vector) - 1
    assert len(results['saccade_y']) == len(fix_vector) - 1
    assert len(results['saccade_lenx']) == len(fix_vector) - 1
    assert len(results['saccade_leny']) == len(fix_vector) - 1
    assert len(results['saccade_theta']) == len(fix_vector) - 1
    assert len(results['saccade_rho']) == len(fix_vector) - 1


def test_cal_vectordifferences(length=5):
    """
    Tests whether vector difference calculation results in matrix of
    expected dimensions.

    Parameters
    -----------
    length: int
       specify the length of the structured array to be used
    """
    data1, data2 = ut.mk_strucarray(length=5)
    Matrix = mp.cal_vectordifferences(data1, data2)
    assert Matrix.shape[0] == length - 1
    assert Matrix.shape[1] == 2 * length - 2


def test_anglesim():
    """
    Tests whether predefined angular relations (quarters of the unit circle)
    yield the expected similarity results.
    """
    M_assignment, path, angles1, angles2, angles3, angles4 = ut.mk_angles()
    a1a2 = mp.cal_angulardifference(angles1, angles2, path, M_assignment)
    a3a4 = mp.cal_angulardifference(angles3, angles4, path, M_assignment)
    assert a1a2 == [0, 1.0459999999999994, 1.5700000000000003, 2.08,
                    3.1400000000000006]
    assert a3a4 == [3.1400000000000006, 2.0971853071795863, 1.5711853071795865,
                    1.0491853071795862, 0.0031853071795859833]


def test_durationsim():
    """
    Tests whether similar but not identical durations result in duration
    similarities greater zero
    """
    M_assignment, path, duration1, duration2 = ut.mk_durs()
    res = mp.cal_durationdifference(duration1, duration2, path, M_assignment)
    assert all(np.asarray(res) >= 0)


def test_compare2matlab():
    """
    compares whether analysis with given inputs yields the same results as a
    calculation with the matlab toolbox
    """
    testfile = os.path.abspath('multimatch/tests/testdata/segment_10_sub-19.tsv')
    testfile2 = os.path.abspath('multimatch/tests/testdata/segment_10_sub-01.tsv')
    data1 = np.recfromcsv(testfile,
                          delimiter='\t',
                          dtype={'names': ('start_x', 'start_y', 'duration'),
                                 'formats': ('f8', 'f8', 'f8')})
    data2 = np.recfromcsv(testfile2,
                          delimiter='\t',
                          dtype={'names': ('start_x', 'start_y', 'duration'),
                                 'formats': ('f8', 'f8', 'f8')})
    res_grouping = mp.docomparison(data1,
                                   data2,
                                   sz=[720, 1280],
                                   grouping=True,
                                   TDir=45.0,
                                   TDur=0.3,
                                   TAmp=146.8)
    res_no_grouping = mp.docomparison(data1,
                                      data2,
                                      sz=[720, 1280],
                                      grouping=False,
                                      TDir=0,
                                      TDur=0,
                                      TAmp=0)
    matlab_grouping = [0.95076, 0.95638, 0.94082, 0.94491, 0.78261]
    matlab_no_grouping = [0.99082, 0.69902, 0.98927, 0.94767, 0.65563]
    from pytest import approx
    assert matlab_grouping == approx(res_grouping[0], abs=1e-5, rel=1e-5)
    assert matlab_no_grouping == approx(res_no_grouping[0], abs=1e-5, rel=1e-5)


def test_too_short_scanpaths():
    """
    if at least one scanpath is too short results should be nan.
    """
    fixvector1 = ut.mk_fix_vector(2)
    fixvector2 = ut.mk_fix_vector(5)
    results = mp.docomparison(fixvector1, fixvector2, sz=[1280, 720], grouping=True, TAmp=8, TDir=8, TDur=8)
    assert ([results[i] == np.nan for i in range(0, len(results))])


# The original use case of this module was for a Masters Thesis that performed
# scanpath comparison on the studyforrest dataset (Hanke et al., 2016). The test
# in this sections utilize (and test) functions that are related to this data
# specifically (found under tests/testdata and in a non-entrypoint script
# multimatch_forrest). Coincidentally, these tests test multimatch functions
# nevertheless, so they'll stay for now..

def test_offsets():
    """
    Tests offset length on example dataset with known results.
    """
    shots = ut.short_shots()
    offsets = ut.create_offsets(shots, dur=5)
    assert len(offsets) == 13

    shots_short = ut.mk_supershort_shots()
    offsets2 = ut.create_offsets(shots_short, dur=2)
    assert len(offsets2) == 0

    shots_long = ut.mk_longershots()
    offsets3 = ut.create_offsets(shots_long, dur=3)
    assert len(offsets3) == len(shots_long)


# disabling those for now -- chunking is not relevant for multimatch anymore

# def test_createchunks(run=1, subj=1):
#     """
#     Tests chunking of studyforrest data into scenes based on shot annotation.
#     Are start and end ids of shots the same length?
#     Does no endid preceed a start id?
#     run: int
#         specify the run (choose between 1 - 8)
#     subj: int
#         specify the subject (example data of lab subject (1)
#         or mri subject (2) available)
#     """
#     data1, data2, shots = ut.same_sample(run, subj)
#     onsets = ut.create_onsets(shots, 3)
#     fixations = ut.preprocess(data1)
#     startid, endid = ut.create_chunks(onsets, fixations, 3)
#     assert len(startid) == len(endid)
#     assert np.all(startid[:-1] <= endid[1:])
#
#
# def test_createoffsetchunks(run=1, subj=1):
#     """
#     Tests chunking of studyforrest data into scenes based on shot annotation.
#     Are start and end ids of shots the same length?
#     Does no endid preceed a start id?
#     run: int
#         specify the run (choose between 1 - 8)
#     subj: int
#         specify the subject (example data of lab subject (1)
#         or mri subject (2) available)
#     """
#     data1, data2, shots = ut.same_sample(run, subj)
#     offsets = ut.create_offsets(shots, 3)
#     fixations = ut.preprocess(data1)
#     startid, endid = ut.create_offsetchunks(offsets, fixations, 3)
#     assert len(startid) == len(endid)
#     assert np.all(startid[:-1] <= endid[1:])


def test_longshot():
    """
    Test whether longshots combines scenes of longer length into the
    correct number of shots.
    """
    shots = ut.short_shots()
    newshots = ut.longshot(shots, group_shots=True, ldur=4.92)
    assert len(newshots) == len(shots) - 2


def test_longshot_nogrouping():
    """
    Tests whether longshots will perform no grouping if dur = None.
    """
    shots = ut.short_shots()
    newshots = ut.longshot(shots, group_shots=False, ldur=0.0)
    assert len(newshots) == len(shots)


def test_closestright():
    """
    Tests whether closestright function works as intended and does not
    unexpectedly return a number left to the reference.
    """
    mylist = np.arange(0, 25)
    for i in range(0, len(mylist)):
        res = ut.takeclosestright(mylist, i)
        if i == 0:
            assert res == 1
        elif i == len(mylist) - 1:
            assert res == 24
        else:
            assert res == i + 1


def test_closestleft():
    """
    Tests whether closestright function works as intended and does not
    unexpectedly return a number right to the reference.
    """
    mylist = np.arange(-1, 25)
    for i in range(0, len(mylist)):
        res = ut.takeclosestleft(mylist, i)
        if i == 0:
            assert res == mylist[0]
        elif i == len(mylist) - 1:
            assert res == 24
        elif i == len(mylist):
            assert res == 25
        else:
            assert res == i - 1
