#!/usr/bin/env python


import numpy as np
import math
import sys
import logging


def cart2pol(x, y):
    """Transform cartesian into polar coordinates.

    :param x: float
    :param y : float

    :return: rho: float, length from (0,0)
    :return: theta: float, angle in radians
    """
    rho = np.sqrt(x ** 2 + y ** 2)
    theta = np.arctan2(y, x)
    return rho, theta


def calcangle(x1, x2):
    """Calculate angle between to vectors (saccades).

    :param: x1, x2: list of float

    :return: angle: float, angle in degrees
    """
    angle = math.degrees(
        math.acos(
            np.dot(x1, x2) / (np.linalg.norm(x1) * np.linalg.norm(x2))))
    return angle


def remodnav_reader(data, screensize, pursuits=False):
    """
    Helper function to read and preprocess REMoDNaV data for use in
    interactive python sessions.

    :param data: path to a REMoDNaV output
    :param screensize: list, screendimensions in x and y direction
    :param pursuits: if True, pursuits will be relabeled to fixations
    """
    from multimatch_gaze.tests import utils as ut
    data = ut.read_remodnav(data)
    # this function can be called without any previous check that
    # screensize are two values, so I'm putting an additional check
    # here
    try:
        assert len(screensize) == 2
    except:
        raise ValueError(
            "Screensize should be the dimensions of the"
            "screen in x and y direction, such as "
            "[1000, 800]. I received {}.".format(
                screensize
            )
        )
    if pursuits:
        data = ut.pursuits_to_fixations(data)
    data = ut.preprocess_remodnav(data, screensize)
    return data


def gen_scanpath_structure(data):
    """Transform a fixation vector into a vector based scanpath representation.

    Takes an nx3 fixation vector (start_x, start_y, duration) in the form of
    of a record array and transforms it into a vector-based scanpath
    representation in the form of a nested dictionary. Saccade starting and
    end points, as well as length in x & y direction, and vector length (theta)
    and direction (rho) are calculated from fixation coordinates as a vector
    representation in 2D space.
    Structure:
    fix --> fixations --> (start_x, start_y, duration)
    sac --> saccades --> (start_x, start_y, lenx, leny, rho, theta)

    :param: data: record array

    :return: eyedata: dict, vector-based scanpath representation
    """

    # everything into a dict
    # keep coordinates and durations of fixations
    fixations = dict(
        x=data['start_x'],
        y=data['start_y'],
        dur=data['duration'],
    )
    # calculate saccade length and angle from vector lengths between fixations
    lenx = np.diff(data['start_x'])
    leny = np.diff(data['start_y'])
    rho, theta = cart2pol(lenx, leny)

    saccades = dict(
        # fixations are the start coordinates for saccades
        x=data[:-1]['start_x'],
        y=data[:-1]['start_y'],
        lenx=lenx,
        leny=leny,
        theta=theta,
        rho=rho,
    )
    return dict(fix=fixations, sac=saccades)


def keepsaccade(i,
                j,
                sim,
                data
                ):
    """
    Helper function for scanpath simplification. If no simplification can be
    performed on a particular saccade, this functions stores the original data.
    :param i: current index
    :param j: current index
    :param sim: dict with current similarities
    :param data: original dict with vector based scanpath representation
    """
    for t, k in (('sac', 'lenx'),
                 ('sac', 'leny'),
                 ('sac', 'x'),
                 ('sac', 'y'),
                 ('sac', 'theta'),
                 ('sac', 'rho'),
                 ('fix', 'dur')):
        sim[t][k].insert(j, data[t][k][i])

    return i + 1, j + 1


def _get_empty_path():
    return dict(
        fix=dict(
            dur=[],
        ),
        sac=dict(
            x=[],
            y=[],
            lenx=[],
            leny=[],
            theta=[],
            # why 'len' here and 'rho' in input data?
            # MIH -> always rho
            #len=[],
            rho=[],
        )
    )


def simlen(path, TAmp, TDur):
    """Simplify scanpaths based on saccadic length.

    Simplify consecutive saccades if their length is smaller than the
    threshold TAmp and the duration of the closest fixations is lower
    than threshold TDur.

    :param: path: dict, output of gen_scanpath_structure
    :param: TAmp: float, length in px
    :param: TDur: float, time in seconds

    :return: eyedata: dict; one iteration of length based simplification
    """
    # shortcuts
    saccades = path['sac']
    fixations = path['fix']

    if len(saccades['x']) < 1:
        return path

    # the scanpath is long enough
    i = 0
    j = 0
    sim = _get_empty_path()
    # while we don't run into index errors
    while i <= len(saccades['x']) - 1:
        # if saccade is the last one
        if i == len(saccades['x']) - 1:
            # and if saccade has a length shorter than the threshold:
            if saccades['rho'][i] < TAmp:
                # and if the fixation duration is short:
                if (fixations['dur'][-1] < TDur) or (fixations['dur'][-2] < TDur):
                    # calculate sum of local vectors for simplification
                    v_x = saccades['lenx'][-2] + saccades['lenx'][-1]
                    v_y = saccades['leny'][-2] + saccades['leny'][-1]
                    rho, theta = cart2pol(v_x, v_y)
                    # save them in the new vectors
                    sim['sac']['lenx'][j - 1] = v_x
                    sim['sac']['leny'][j - 1] = v_y
                    sim['sac']['theta'][j - 1] = theta
                    sim['sac']['rho'][j - 1] = rho
                    sim['fix']['dur'].insert(j, fixations['dur'][i - 1])
                    j -= 1
                    i += 1
                # if fixation duration is longer than the threshold:
                else:
                    # insert original event data in new list -- no
                    # simplification
                    i, j = keepsaccade(i, j, sim, path)
            # if saccade does NOT have a length shorter than the threshold:
            else:
                # insert original path in new list -- no simplification
                i, j = keepsaccade(i, j, sim, path)
        # if saccade is not the last one
        else:
            # and if saccade has a length shorter than the threshold
            if (saccades['rho'][i] < TAmp) and (i < len(saccades['x']) - 1):
                # and if fixation durations are short
                if (fixations['dur'][i + 1] < TDur) or \
                        (fixations['dur'][i] < TDur):
                    # calculate sum of local vectors in x and y length for
                    # simplification
                    v_x = saccades['lenx'][i] + saccades['lenx'][i + 1]
                    v_y = saccades['leny'][i] + saccades['leny'][i + 1]
                    rho, theta = cart2pol(v_x, v_y)
                    # save them in the new vectors
                    sim['sac']['lenx'].insert(j, v_x)
                    sim['sac']['leny'].insert(j, v_y)
                    sim['sac']['x'].insert(j, saccades['x'][i])
                    sim['sac']['y'].insert(j, saccades['y'][i])
                    sim['sac']['theta'].insert(j, theta)
                    sim['sac']['rho'].insert(j, rho)
                    # add the old fixation duration
                    sim['fix']['dur'].insert(j, fixations['dur'][i])
                    i += 2
                    j += 1
                # if fixation durations longer than the threshold
                else:
                    # insert original path in new lists -- no simplification
                    i, j = keepsaccade(i, j, sim, path)
            # if saccade does NOT have a length shorter than the threshold:
            else:
                # insert original path in new list -- no simplification
                i, j = keepsaccade(i, j, sim, path)
    # append the last fixation duration
    sim['fix']['dur'].append(fixations['dur'][-1])

    return sim


def simdir(path,
           TDir,
           TDur
           ):
    """Simplify scanpaths based on angular relations between saccades (direction).

    Simplify consecutive saccades if the angle between them is smaller than the
    threshold TDir and the duration of the intermediate fixations is lower
    than threshold TDur.

    :param: path: dict, output of gen_scanpath_structure
    :param: TDir: float, angle in degrees
    :param: TDur: float, time in seconds

    :return: eyedata: dict, one iteration of direction based simplification
    """
    # shortcuts
    saccades = path['sac']
    fixations = path['fix']

    if len(saccades['x']) < 1:
        return path
    # the scanpath is long enough
    i = 0
    j = 0
    sim = _get_empty_path()
    # while we don't run into index errors
    while i <= len(saccades['x']) - 1:
        if i < len(saccades['x']) - 1:
            # lets check angles
            v1 = [saccades['lenx'][i], saccades['leny'][i]]
            v2 = [saccades['lenx'][i + 1], saccades['leny'][i + 1]]
            angle = calcangle(v1, v2)
        else:
            # an angle of infinite size won't go into any further loop
            angle = float('inf')
        # if the angle is smaller than the threshold and its not the last saccade
        if (angle < TDir) & (i < len(saccades['x']) - 1):
            # if the fixation duration is short:
            if fixations['dur'][i + 1] < TDur:
                # calculate the sum of local vectors
                v_x = saccades['lenx'][i] + saccades['lenx'][i + 1]
                v_y = saccades['leny'][i] + saccades['leny'][i + 1]
                rho, theta = cart2pol(v_x, v_y)
                # save them in the new vectors
                sim['sac']['lenx'].insert(j, v_x)
                sim['sac']['leny'].insert(j, v_y)
                sim['sac']['x'].insert(j, saccades['x'][i])
                sim['sac']['y'].insert(j, saccades['y'][i])
                sim['sac']['theta'].insert(j, theta)
                sim['sac']['rho'].insert(j, rho)
                # add the fixation duration
                sim['fix']['dur'].insert(j, fixations['dur'][i])
                i += 2
                j += 1
            else:
                # insert original data in new list -- no simplification
                i, j = keepsaccade(i, j, sim, path)
        # elif the angle is smaller than the threshold, but its the LAST saccade:
        ## Testing revealed that we never actually get here -- because for the
        ## last saccade, the angle is inf. This however, is how it seems to be
        ## implemented in the original toolbox.
        ##  TODO: ponder whether to keep exact original (dys)functionality here
        # elif (angle < TDir) & (i == len(saccades['x']) - 1):
        #     print("step 1", angle, i)
        #     # if the fixation duration is short:
        #     if fixations['dur'][i + 1] < TDur:
        #         # calculate sum of local vectors
        #         print("TRIGGERED")
        #         v_x = saccades['lenx'][i - 2] + saccades['lenx'][i - 1]
        #         v_y = saccades['leny'][i - 2] + saccades['leny'][i - 1]
        #         rho, theta = cart2pol(v_x, v_y)
        #         # save them in new vectors
        #         sim['sac']['lenx'][j - 1] = v_x
        #         sim['sac']['leny'][j - 1] = v_y
        #         sim['sac']['theta'][j - 1] = theta
        #         sim['sac']['len'][j - 1] = rho
        #         sim['fix']['dur'].insert(j, fixations['dur'][-1] + (fixations['dur'][i] / 2))
        #         j -= 1
        #         i += 1
        #     # if fixation duration is longer than the threshold:
        #     else:
        #         # insert original path in new list -- no simplification
        #         i, j = keepsaccade(i, j, sim, path)
        # else (the angle is larger than the threshold)
        else:
            # insert original path in new list -- no simplification
            i, j = keepsaccade(i, j, sim, path)
    # now append the last fixation duration
    sim['fix']['dur'].append(fixations['dur'][-1])

    return sim


def simplify_scanpath(path,
                      TAmp,
                      TDir,
                      TDur
                      ):
    """Simplify scanpaths until no further simplification is possible.

    Loops over simplification functions simdir and simlen until no
    further simplification of the scanpath is possible.

    :param: path: dict, vector based scanpath representation,
                  output of gen_scanpath_structure
    :param: TAmp: float, length in px
    :param: TDir: float, angle in degrees
    :param: TDur: float, duration in seconds

    :return: eyedata: dict, simplified vector-based scanpath representation
    """
    looptime = 0
    while True:
        path = simdir(path, TDir, TDur)
        path = simlen(path, TAmp, TDur)
        looptime += 1
        if looptime == len(path['fix']['dur']):
            return path


def cal_vectordifferences(path1,
                          path2
                          ):
    """Create matrix of vector-length differences of all vector pairs

    Create M, a Matrix with all possible saccade-length differences between
    saccade pairs.

    :param: path1, path2: dicts, vector-based scanpath representations

    :return: M: array-like
        Matrix of vector length differences

    """
    # take length in x and y direction of both scanpaths
    x1 = np.asarray(path1['sac']['lenx'])
    x2 = np.asarray(path2['sac']['lenx'])
    y1 = np.asarray(path1['sac']['leny'])
    y2 = np.asarray(path2['sac']['leny'])
    # initialize empty list for rows, will become matrix to store sacc-length
    # pairings
    rows = []
    # calculate saccade length differences, vectorized
    for i in range(0, len(x1)):
        x_diff = abs(x1[i] * np.ones(len(x2)) - x2)
        y_diff = abs(y1[i] * np.ones(len(y2)) - y2)
        # calc final length from x and y lengths, append, stack into matrix M
        rows.append(np.asarray(np.sqrt(x_diff ** 2 + y_diff ** 2)))
    M = np.vstack(rows)
    return M


def createdirectedgraph(scanpath_dim,
                        M,
                        M_assignment
                        ):
    """Create a directed graph:
    The data structure of the result is a nested dictionary such as
    weightedGraph = {0 : {1:259.55, 15:48.19, 16:351.95},
    1 : {2:249.354, 16:351.951, 17:108.97},
    2 : {3:553.30, 17:108.97, 18:341.78}, ...}

    It defines the possible nodes to reach from a particular node, and the weight that
    is associated with the path to each of the possible nodes.

    :param: scanpath_dim: list, shape of matrix M
    :param: M: array-like, matrix of vector length differences
    :param: M_assignment: array-like, Matrix, arranged with values from 0 to number of entries in M

    :return: weighted graph: dict, Dictionary within a dictionary pairing weights (distances) with
            node-pairings

    """

    # initialize dictionary for neighbouring vertices and edge weights
    adjacent = {}
    weight = {}
    # loop through every node rowwise
    for i in range(0, scanpath_dim[0]):
        # loop through every node columnwise
        for j in range(0, scanpath_dim[1]):
            currentNode = i * scanpath_dim[1] + j
            # if in the last (bottom) row, only go right
            if (i == scanpath_dim[0] - 1) & (j < scanpath_dim[1] - 1):
                adjacent[M_assignment[i, j]] = [currentNode + 1]
                weight[M_assignment[i, j]] = [M[i, j + 1]]
            # if in the last (rightmost) column, only go down
            elif (i < scanpath_dim[0] - 1) & (j == scanpath_dim[1] - 1):
                adjacent[M_assignment[i, j]] = [currentNode + scanpath_dim[1]]
                weight[M_assignment[i, j]] = [M[i + 1, j]]
            # if in the last (bottom-right) vertex, do not move any further
            elif (i == scanpath_dim[0] - 1) & (j == scanpath_dim[1] - 1):
                adjacent[M_assignment[i, j]] = [currentNode]
                weight[M_assignment[i, j]] = [0]
            # anywhere else, move right, down and down-right.
            else:
                adjacent[M_assignment[i, j]] = [currentNode + 1,
                                                currentNode + scanpath_dim[1],
                                                currentNode + scanpath_dim[1] + 1]
                weight[M_assignment[i, j]] = [M[i, j + 1],
                                              M[i + 1, j],
                                              M[i + 1, j + 1]]
    # create ascending list ranging from first to last node - this
    #  will be the first key in the nested dict
    Startnodes = range(0, scanpath_dim[0] * scanpath_dim[1])
    # initialize list with adjacent nodes (adjacent to each startnode)
    # and the weights associated with the paths between them
    weightedEdges = [
        dict(zip(a, w)) for a, w in zip(adjacent.values(), weight.values())
    ]
    # initialize final dictionary
    weightedGraph = dict(zip(Startnodes, weightedEdges))
    return weightedGraph


def dijkstra(weightedGraph,
             start,
             end
             ):
    """Implementation of Dijkstra algorithm:
    Use the dijkstra algorithm to find the shortest path through a directed
    graph (weightedGraph) from start to end.

    :param: weightedGraph: dict, dictionary within a dictionary pairing weights (distances) with
            node-pairings
    :param: start: int, starting point of path, should be 0
    :param: end: int, end point of path, should be (n, m) of Matrix M

    :return: path: array, indices of the shortest path, i.e. best-fitting saccade pairs
    :return: dist: float, sum of weights

    """

    # initialize empty dictionary to hold distances
    dist = {}
    # inialize list of vertices in the path to current vertex (predecessors)
    pred = {}
    # where do I need to go still?
    to_assess = weightedGraph.keys()
    for node in weightedGraph:
        # set inital distances to infinity
        dist[node] = float('inf')
        # no node has any predecessors yet
        pred[node] = None
    # initialize list to be filled with final distances(weights) of nodes
    sp_set = []
    # the starting node gets a weight of 0 to make sure to start there
    dist[start] = 0
    # continue the algorithm as long as there are still unexplored nodes
    while len(sp_set) < len(to_assess):
        still_in = {node: dist[node] for node in [node for node in to_assess if
                                                  node not in sp_set]}
        # find adjacent node with minimal weight and append to sp_set
        closest = min(still_in, key=dist.get)
        sp_set.append(closest)
        for node in weightedGraph[closest]:
            if dist[node] > dist[closest] + weightedGraph[closest][node]:
                dist[node] = dist[closest] + weightedGraph[closest][node]
                pred[node] = closest
    # append endnode to list path
    path = [end]
    # append contents of pred in reversed order to path
    while start not in path:
        path.append(pred[path[-1]])
    # return path in reverse order (begin to end) and final distance
    return path[::-1], dist[end]


def cal_angulardifference(data1,
                          data2,
                          path,
                          M_assignment
                          ):
    """Calculate angular similarity of two scanpaths:

    :param: data1: dict; contains vector-based scanpath representation of the
        first scanpath
    :param: data2: dict, contains vector-based scanpath representation of the
        second scanpath
    :param: path: array,
        indices for the best-fitting saccade pairings between scanpaths
    :param: M_assignment: array-like, Matrix arranged with values from 0 to number of entries in
        M, the matrix of vector length similarities

    :return: anglediff: array of floats, angular differences between pairs of saccades
        of two scanpaths

    """
    # get the angle between saccades from the scanpaths
    theta1 = data1['sac']['theta']
    theta2 = data2['sac']['theta']
    # initialize list to hold individual angle differences
    anglediff = []
    # calculate angular differences between the saccades along specified path
    for p in path:
        # which saccade indices correspond to path?
        i, j = np.where(M_assignment == p)
        # extract the angle
        spT = [theta1[i.item()], theta2[j.item()]]
        for t in range(0, len(spT)):
            # get results in range -pi, pi
            if spT[t] < 0:
                spT[t] = math.pi + (math.pi + spT[t])
        spT = abs(spT[0] - spT[1])
        if spT > math.pi:
            spT = 2 * math.pi - spT
        anglediff.append(spT)
    return anglediff


def cal_durationdifference(data1,
                           data2,
                           path,
                           M_assignment
                           ):
    """Calculate similarity of two scanpaths fixation durations.

    :param: data1: array-like
        dict, contains vector-based scanpath representation of the
        first scanpath
    :param: data2: array-like
        dict, contains vector-based scanpath representation of the
        second scanpath
    :param: path: array
        indices for the best-fitting saccade pairings between scanpaths
    :param: M_assignment: array-like
         Matrix, arranged with values from 0 to number of entries in M, the
         matrix of vector length similarities

    :return: durdiff: array of floats,
        array of fixation duration differences between pairs of saccades from
        two scanpaths

    """
    # get the duration of fixations in the scanpath
    dur1 = data1['fix']['dur']
    dur2 = data2['fix']['dur']
    # initialize list to hold individual duration differences
    durdiff = []
    # calculation fixation duration differences between saccades along path
    for p in path:
        # which saccade indices correspond to path?
        i, j = np.where(M_assignment == p)
        maxlist = [dur1[i.item()], dur2[j.item()]]
        # compute abs. duration diff, normalize by largest duration in pair
        durdiff.append(abs(dur1[i.item()] -
                           dur2[j.item()]) / abs(max(maxlist)))
    return durdiff


def cal_lengthdifference(data1,
                         data2,
                         path,
                         M_assignment
                         ):
    """Calculate length similarity of two scanpaths.

    :param: data1: array-like
        dict, contains vector-based scanpath representation of the
        first scanpath
    :param: data2: array-like
        dict, contains vector-based scanpath representation of the
        second scanpath
    :param: path: array
        indices for the best-fitting saccade pairings between scanpaths
    :param: M_assignment: array-like
         Matrix, arranged with values from 0 to number of entries in M, the
         matrix of vector length similarities

    :return: lendiff: array of floats
        array of length difference between pairs of saccades of two scanpaths

    """
    # get the saccade lengths rho
    len1 = np.asarray(data1['sac']['rho'])
    len2 = np.asarray(data2['sac']['rho'])
    # initialize list to hold individual length differences
    lendiff = []
    # calculate length differences between saccades along path
    for p in path:
        i, j = np.where(M_assignment == p)
        lendiff.append(abs(len1[i] - len2[j]))
    return lendiff


def cal_positiondifference(data1,
                           data2,
                           path,
                           M_assignment
                           ):
    """Calculate position similarity of two scanpaths.

    :param: data1: array-like
        dict, contains vector-based scanpath representation of the
        first scanpath
    :param: data2: array-like
        dict, contains vector-based scanpath representation of the
        second scanpath
    :param: path: array
        indices for the best-fitting saccade pairings between scanpaths
    :param: M_assignment: array-like
         Matrix, arranged with values from 0 to number of entries in M, the
         matrix of vector length similarities

    :return: posdiff: array of floats
        array of position differences between pairs of saccades
        of two scanpaths

    """
    # get the x and y coordinates of points between saccades
    x1 = np.asarray(data1['sac']['x'])
    x2 = np.asarray(data2['sac']['x'])
    y1 = np.asarray(data1['sac']['y'])
    y2 = np.asarray(data2['sac']['y'])
    # initialize list to hold individual position differences
    posdiff = []
    # calculate position differences along path
    for p in path:
        i, j = np.where(M_assignment == p)
        posdiff.append(math.sqrt((x1[i.item()] - x2[j.item()]) ** 2 +
                                 (y1[i.item()] - y2[j.item()]) ** 2))
    return posdiff


def cal_vectordifferencealongpath(data1,
                                  data2,
                                  path,
                                  M_assignment
                                  ):
    """Calculate vector similarity of two scanpaths.

    :param: data1: array-like
        dict, contains vector-based scanpath representation of the
        first scanpath
    :param: data2: array-like
        dict, contains vector-based scanpath representation of the
        second scanpath
    :param: path: array-like
        array of indices for the best-fitting saccade pairings between scan-
        paths
    :param: M_assignment: array-like
         Matrix, arranged with values from 0 to number of entries in M, the
         matrix of vector length similarities

    :return: vectordiff: array of floats
            array of vector differences between pairs of saccades of two scanpaths

    """
    # get the saccade lengths in x and y direction of both scanpaths
    x1 = np.asarray(data1['sac']['lenx'])
    x2 = np.asarray(data2['sac']['lenx'])
    y1 = np.asarray(data1['sac']['leny'])
    y2 = np.asarray(data2['sac']['leny'])
    # initialize list to hold individual vector differences
    vectordiff = []
    # calculate vector differences along path
    # TODO look at this again, should be possible simpler
    for p in path:
        i, j = np.where(M_assignment == p)
        vectordiff.append(np.sqrt((x1[i.item()] - x2[j.item()]) ** 2 +
                                  (y1[i.item()] - y2[j.item()]) ** 2))
    return vectordiff


def getunnormalised(data1,
                    data2,
                    path,
                    M_assignment
                    ):
    """Calculate unnormalised similarity measures.

    Calls the five functions to create unnormalised similarity measures for
    each of the five similarity dimensions. Takes the median of the resulting
    similarity values per array.

    :param: data1: array-like
        dict, contains vector-based scanpath representation of the
        first scanpath
    :param: data2: array-like
        dict, contains vector-based scanpath representation of the
        second scanpath
    :param: path: array
        indices for the best-fitting saccade pairings between scanpaths
    :param: M_assignment: array-like
         Matrix, arranged with values from 0 to number of entries in M, the
         matrix of vector length similarities

    :return: unnormalised: array
        array of unnormalised similarity measures on five dimensions

    >>> unorm_res = getunnormalised(scanpath_rep1, scanpath_rep2, path, M_assignment)
    """
    return [
        np.median(fx(data1, data2, path, M_assignment))
        for fx in (cal_vectordifferencealongpath,
                   cal_angulardifference,
                   cal_lengthdifference,
                   cal_positiondifference,
                   cal_durationdifference)
    ]


def normaliseresults(unnormalised, screensize):
    """Normalize similarity measures.

    Vector similarity is normalised against two times screen diagonal,
    the maximum theoretical distance.
    Direction similarity is normalised against pi.
    Length Similarity is normalised against screen diagonal.
    Position Similarity and Duration Similarity are already normalised.

    :param: unnormalised: array
        array of unnormalised similarity measures,
        output of getunnormalised()

    :return: normalresults: array
        array of normalised similarity measures

    >>> normal_res = normaliseresults(unnormalised, screensize)
    """
    # normalize vector similarity against two times screen diagonal, the maximum
    # theoretical distance
    VectorSimilarity = 1 - unnormalised[0] / (2 * math.sqrt(screensize[0] ** 2 + screensize[1] ** 2))
    # normalize against pi
    DirectionSimilarity = 1 - unnormalised[1] / math.pi
    # normalize against screen diagonal
    LengthSimilarity = 1 - unnormalised[2] / math.sqrt(screensize[0] ** 2 + screensize[1] ** 2)
    PositionSimilarity = 1 - unnormalised[3] / math.sqrt(screensize[0] ** 2 + screensize[1] ** 2)
    # no normalisazion necessary, already done
    DurationSimilarity = 1 - unnormalised[4]
    normalresults = [VectorSimilarity, DirectionSimilarity, LengthSimilarity,
                     PositionSimilarity, DurationSimilarity]
    return normalresults


def docomparison(fixation_vectors1,
                 fixation_vectors2,
                 screensize,
                 grouping=False,
                 TDir=0.0,
                 TDur=0.0,
                 TAmp=0.0
                 ):
    """Compare two scanpaths on five similarity dimensions.


    :param: fixation_vectors1: array-like n x 3 fixation vector of one scanpath
    :param: fixation_vectors2: array-like n x 3 fixation vector of one scanpath
    :param: screensize: list, screen dimensions in px.
    :param: grouping: boolean, if True, simplification is performed based on thresholds TAmp,
        TDir, and TDur. Default: False
    :param: TDir: float, Direction threshold, angle in degrees. Default: 0.0
    :param: TDur: float,  Duration threshold, duration in seconds. Default: 0.0
    :param: TAmp: float, Amplitude threshold, length in px. Default: 0.0

    :return: scanpathcomparisons: array
        array of 5 scanpath similarity measures. Vector (Shape), Direction
        (Angle), Length, Position, and Duration. 1 means absolute similarity, 0 means
        lowest similarity possible.

    >>> results = docomparison(fix_1, fix_2, screensize = [1280, 720], grouping = True, TDir = 45.0, TDur = 0.05, TAmp = 150)
    >>> print(results)
    >>> [[0.95075847681364678, 0.95637548674423822, 0.94082367355291008, 0.94491164030498609, 0.78260869565217384]]
    """
    # check if fixation vectors/scanpaths are long enough
    if (len(fixation_vectors1) >= 3) & (len(fixation_vectors2) >= 3):
        # get the data into a geometric representation
        path1 = gen_scanpath_structure(fixation_vectors1)
        path2 = gen_scanpath_structure(fixation_vectors2)
        if grouping:
            # simplify the data
            path1 = simplify_scanpath(path1, TAmp, TDir, TDur)
            path2 = simplify_scanpath(path2, TAmp, TDir, TDur)
        # create M, a matrix of all vector pairings length differences (weights)
        M = cal_vectordifferences(path1, path2)
        # initialize a matrix of size M for a matrix of nodes
        scanpath_dim = np.shape(M)
        M_assignment = np.arange(scanpath_dim[0] * scanpath_dim[1]).reshape(scanpath_dim[0], scanpath_dim[1])
        # create a weighted graph of all possible connections per Node, and their weight
        weightedGraph = createdirectedgraph(scanpath_dim, M, M_assignment)
        # find the shortest path (= lowest sum of weights) through the graph
        path, dist = dijkstra(weightedGraph, 0, scanpath_dim[0] * scanpath_dim[1] - 1)
        # compute similarities on alinged scanpaths and normalize them
        unnormalised = getunnormalised(path1, path2, path, M_assignment)
        normal = normaliseresults(unnormalised, screensize)
        return normal
    # return nan as result if at least one scanpath it too short
    else:
        return np.repeat(np.nan, 5)


def parse_args(args):
    """Argument parse for command line invocation

    Turned it into a function to make testing easier.

    :param args: [command line] arguments
    :return: argument parser
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog='multimatch_gaze',
        description='{}'.format(
            main.__doc__
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        'input1', metavar='<datafile>',
        help="""Fixation data of scanpath 1. Should be a tab separated
        file with columns corresponding to x-coordinates ('start_x'),
        y-coordinates ('start_y'), and fixation duration ('duration')
        in seconds.""")
    parser.add_argument(
        'input2', metavar='<datafile>',
        help="""Fixation data of scanpath 2. Should be a tab separated
        file with columns corresponding to x-coordinates ('start_x'),
        y-coordinates ('start_y'), and fixation duration ('duration')
        in seconds.""")
    parser.add_argument(
        'screensize',  metavar='<screensize>',
        nargs='+',
        help="""screensize: Resolution of screen in px, should be supplied as
        1000 800 for a screen of resolution [1000, 800]. This parameter is
        necessary to correctly normalize Length, Position, and Vector similarity
        to range [0, 1].""")
    parser.add_argument(
        '--direction-threshold', type=float, metavar='<TDir>', default=0.0,
        help="""Threshold for direction based grouping in degree (example: 45.0).
        Two consecutive saccades with an angle below TDir and short fixations will
        be grouped together to reduce scanpath complexity. If 0: no
        simplification will be performed.
        Default: 0 (no simplification)""")
    parser.add_argument(
        '--amplitude-threshold', type=float, metavar='<TAmp>', default=0.0,
        help="""Threshold for amplitude based grouping in pixel (example: 140.0).
        Two consecutive saccades shorter than TAmp and short fixations will be
        grouped together to reduce scanpath complexity. If 0: no simplification
        will be performed.
        Default: 0 (no simplification)""")
    parser.add_argument(
        '--duration-threshold', type=float, metavar='<TDur>', default=0.0,
        help="""Threshold for fixation duration during amplitude and direction
        based grouping, in seconds (example: 0.1).
        Default: 0 (no simplification)""")
    parser.add_argument(
        '-o', '--output-type',
        help="""Specify output format of the results: "hr", "single-row"
        or "single-del".
        <hr>: the most Human Readable option, will print dimension
        and value row-wise to the terminal.
        <single-row>: useful to collate results in a table, will print the
        values in a tab-seperated, single string.
        <single-del>: print dimension and value separated with a single
        delimiter (tab), row-wise, without whitespace. Useful to pick a selection
        of scores, split by a single tab, without worrying about whitespace
        default: hr""",
        default = 'hr')
    parser.add_argument(
        '--remodnav', default=False, action='store_true',
        help="""If the input files are output of the REMoDNaV algorithm, and
        the --remodnav parameter is given, multimatch-gaze will read in the
        REMoDNaV data natively. default: False""")
    parser.add_argument(
        '--pursuit', choices=('discard', 'keep'),
        help="""IF the --remodnav parameter is given: Which action to take to
        deal with results? Chose from: 'discard', 'keep'.
        Discard will discard any pursuit event.
        Keep will keep start and end points of pursuits in the
        gaze path.""")

    return parser.parse_args(args)


def main(args=None):
    """Multimatch-gaze: Scanpath comparison in Python.

     Multimatch-gaze is a Python-based reimplementation of the MultiMatch method
     for scanpath comparison (Jarodzka et al., 2010; Dewhurst et al., 2012).
     Based on A) two tab-separated scanpath input files that contain the start x-
     and y-coordinates of fixations and their durations, and B) the screensize in
     pixel, multimatch_gaze calculates the similarity of the provided scanpaths
     on the five dimensions 'shape', 'direction', 'fixation duration', 'length',
     and position (normed to range [0, 1]).
     Scanpath simplification based on angular relation or length is possible on demand.

     For further information, please see https://multimatch_gaze.readthedocs.io/en/latest/.


    """
    # I want to give infos to the user in the command line, but it shouldn't
    # go to stdout -- that would make collation in a table horrible.
    logging.basicConfig(
        format='%(levelname)s:%(message)s',
        level=logging.INFO)
    # I'm sure this function parameter is ugly -- I'm trying to test main with
    # my unit test, in which I need to pass the args...
    if not args:
        args = parse_args(sys.argv[1:])

    screensize = [float(i) for i in args.screensize]
    if len(screensize) != 2:
        raise ValueError(
            'I expected two floats after for the positional'
            'screensize argument, such as 1280 720. '
            'However, I got {}. Please provide the screensize'
            'in pixel')

    if args.remodnav:
        from multimatch_gaze.tests import utils as ut
        # read in the remodnav data
        data1 = ut.read_remodnav(args.input1)
        data2 = ut.read_remodnav(args.input2)

        if args.pursuit == 'keep':
            data1 = ut.pursuits_to_fixations(data1)
            data2 = ut.pursuits_to_fixations(data2)
            #print("Triggered")
            #import pdb; pdb.set_trace()

        data1 = ut.preprocess_remodnav(data1, screensize)
        data2 = ut.preprocess_remodnav(data2, screensize)
    else:
        data1 = np.recfromcsv(args.input1,
                              delimiter='\t',
                              dtype={'names': ('start_x', 'start_y', 'duration'),
                                     'formats': ('f8', 'f8', 'f8')},
                              usecols=(0, 1, 2)
                              )
        data2 = np.recfromcsv(args.input2,
                              delimiter='\t',
                              dtype={'names': ('start_x', 'start_y', 'duration'),
                                     'formats': ('f8', 'f8', 'f8')},
                              usecols=(0, 1, 2)
                              )

    TDir = args.direction_threshold
    TAmp = args.amplitude_threshold
    TDur = args.duration_threshold

    if (TDir != 0) and (TAmp != 0):
        grouping = True
        # give information about the specified analysis, but to stderr
        logging.info(
            'Scanpath comparison is done with simplification. Two consecutive '
            'saccades shorter than {}px and '
            'with an angle smaller than {} degrees are grouped together if '
            'intermediate fixations are shorter '
            'than {} seconds.'.format(TAmp, TDir, TDur))
    else:
        grouping = False
        logging.info(
            'Scanpath comparison is done without any simplification.')

    allowed_output = ['hr', 'single-row', 'single-del']
    output = args.output_type if args.output_type in allowed_output else False

    if not output:
        raise ValueError(
                "I expected an output type specification of 'hr', 'single-row'"
                " or 'single-del', supplied as a string (as in -o 'single-row')."
                " However, I got '{}' instead.".format(args.output_type)
                )
    result = docomparison(data1,
                          data2,
                          screensize=screensize,
                          grouping=grouping,
                          TDir=TDir,
                          TDur=TDur,
                          TAmp=TAmp)

    for i, label in enumerate(('Vector',
                               'Direction',
                               'Length',
                               'Position',
                               'Duration')):
        if output == 'hr':
            print('{} similarity = {}'.format(label, result[i]))
        elif output == 'single-del':
            print('{}\t{}\t'.format(label, result[i]))

    if output == 'single-row':
        print('{}\t{}\t{}\t{}\t{}\t'.format(result[0],
                                            result[1],
                                            result[2],
                                            result[3],
                                            result[4]))
if __name__ == '__main__':

    # execution
    main()
