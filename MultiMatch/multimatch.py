#!/usr/bin/python3


import numpy as np
import math


def cart2pol(x, y):
    """Transform cartesian into polar coordinates.

    Parameters
    ---------
    x, y : float

    Returns
    ---------
    rho: float
        length from (0,0)
    theta: float
        angle in radians
    """
    rho = np.sqrt(x ** 2 + y ** 2)
    theta = np.arctan2(y, x)
    return rho, theta


def calcangle(x1, x2):
    """Calculate angle between to vectors (saccades).

    Parameters
    ---------
    x1, x2: list of float

    Returns
    ---------
    angle: float
        angle in degrees
    """
    angle = math.degrees(
        math.acos(
            np.dot(x1, x2) / (np.linalg.norm(x1) * np.linalg.norm(x2))))
    return angle


def gen_scanpath_structure(data):
    """Transform a fixation vector into a vector based scanpath representation.

    Takes an nx3 fixation vector (start_x, start_y, duration) in the form of
    of a record array and transforms it into appropriate vectorbased scanpath
    representation. Indices are as follows:
    0: fixation_x
    1: fixation_y
    2: fixation_dur
    3: saccade_x
    4: saccade_y
    5: saccade_lenx
    6: saccade_leny
    7: saccade_theta
    8: saccade_rho

    Parameters
    ---------
    data: record array

    Returns
    ---------
    eyedata: array-like
        list of lists, vector-based scanpath representation
    """

    # initialize empty lists
    fixation_x = []
    fixation_y = []
    fixation_dur = []
    saccade_x = []
    saccade_y = []
    saccade_lenx = []
    saccade_leny = []
    saccade_theta = []
    saccade_rho = []
    # get the number of rows
    length = np.shape(data)[0]
    # keep coordinates and durations of fixations
    for i in range(0, length):
        fixation_x.append(data[i]['start_x'])
        fixation_y.append(data[i]['start_y'])
        fixation_dur.append(data[i]['duration'])
    # fixations as start coordinates for saccades (ignores PSOs and the like)
    for i in range(0, length - 1):
        saccade_x.append(data[i]['start_x'])
        saccade_y.append(data[i]['start_y'])
    # calculate saccade length and angle
    for i in range(1, length):
        saccade_lenx.append(fixation_x[i] - saccade_x[i - 1])
        saccade_leny.append(fixation_y[i] - saccade_y[i - 1])
        rho, theta = cart2pol(saccade_lenx[i - 1], saccade_leny[i - 1])
        saccade_rho.append(rho)
        saccade_theta.append(theta)
    # append everything. Eyedata is a list of lists.
    eyedata = [
               fixation_x,
               fixation_y,
               fixation_dur,
               saccade_x,
               saccade_y,
               saccade_lenx,
               saccade_leny,
               saccade_theta,
               saccade_rho
               ]
    return eyedata


def simlen(data, TAmp, TDur):
    """Simplify scanpaths based on saccadic length.

    Simplify consecutive saccades if their length is smaller than the
    threshold TAmp and the duration of the closest fixations is lower
    than threshold TDur.

    Parameters
    ---------
    data: array-like
        list of lists, output of gen_scanpath_structure
    TAmp: float
        length in px
    TDur: float
        time in seconds

    Returns
    ---------
    eyedata: list of lists, one iteration of length based simplification
    """

    if len(data[3]) < 1:
        return data
    # if the scanpath is long enough
    else:
        i = 0
        j = 0
        # initialize empty lists
        sim_dur = []
        sim_x = []
        sim_y = []
        sim_lenx = []
        sim_leny = []
        sim_theta = []
        sim_len = []
        # while we don't run into index errors
        while i <= len(data[3]) - 1:
            # if saccade is the last saccade
            if i == len(data[3]) - 1:
                # if saccade has short length:
                if data[8][i] < TAmp:
                    # if fixation duration is short:
                    if (data[2][-1] < TDur) or (data[2][-2] < TDur):
                        # calculate sum of local vectors
                        v_x = data[5][-2] + data[5][-1]
                        v_y = data[6][-2] + data[6][-1]
                        rho, theta = cart2pol(v_x, v_y)
                        # save them in the new vectors
                        sim_lenx[j - 1] = v_x
                        sim_leny[j - 1] = v_y
                        sim_theta[j - 1] = theta
                        sim_len[j - 1] = rho
                        sim_dur.insert(j, data[2][i - 1])
                        j -= 1
                        i += 1
                    # if fixation duration is long:
                    else:
                        # insert original data in new list
                        sim_lenx.insert(j, data[5][i])
                        sim_leny.insert(j, data[6][i])
                        sim_x.insert(j, data[3][i])
                        sim_y.insert(j, data[4][i])
                        sim_theta.insert(j, data[7][i])
                        sim_len.insert(j, data[8][i])
                        sim_dur.insert(j, data[2][i])
                        i += 1
                        j += 1
                # if saccade doesn't have short length:
                else:
                    # insert original data in new list
                    sim_lenx.insert(j, data[5][i])
                    sim_leny.insert(j, data[6][i])
                    sim_x.insert(j, data[3][i])
                    sim_y.insert(j, data[4][i])
                    sim_theta.insert(j, data[7][i])
                    sim_len.insert(j, data[8][i])
                    sim_dur.insert(j, data[2][i])
                    i += 1
                    j += 1
            # if saccade is not the last one
            else:
                # if saccade has short length
                if (data[8][i] < TAmp) and (i < len(data[3]) - 1):
                    # if fixation durations are short
                    if (data[2][i + 1] < TDur) or (data[2][i] < TDur):
                        # calculate sum of local vectors in x and y length
                        v_x = data[5][i] + data[5][i + 1]
                        v_y = data[6][i] + data[6][i + 1]
                        rho, theta = cart2pol(v_x, v_y)
                        # save them in the new vectors
                        sim_lenx.insert(j, v_x)
                        sim_leny.insert(j, v_y)
                        sim_x.insert(j, data[3][i])
                        sim_y.insert(j, data[4][i])
                        sim_theta.insert(j, theta)
                        sim_len.insert(j, rho)
                        # add the old fixation duration
                        sim_dur.insert(j, data[2][i])
                        i += 2
                        j += 1
                    # if fixation durations are long
                    else:
                        # insert original data in new lists
                        sim_lenx.insert(j, data[5][i])
                        sim_leny.insert(j, data[6][i])
                        sim_x.insert(j, data[3][i])
                        sim_y.insert(j, data[4][i])
                        sim_theta.insert(j, data[7][i])
                        sim_len.insert(j, data[8][i])
                        sim_dur.insert(j, data[2][i])
                        j += 1
                        i += 1
                # if saccade doesn't have short length
                else:
                    # insert original data in new list
                    sim_lenx.insert(j, data[5][i])
                    sim_leny.insert(j, data[6][i])
                    sim_x.insert(j, data[3][i])
                    sim_y.insert(j, data[4][i])
                    sim_theta.insert(j, data[7][i])
                    sim_len.insert(j, data[8][i])
                    sim_dur.insert(j, data[2][i])
                    i += 1
                    j += 1
    sim_dur.append(data[2][-1])
    eyedata = [
               [],
               [],
               sim_dur,
               sim_x,
               sim_y,
               sim_lenx,
               sim_leny,
               sim_theta,
               sim_len
               ]
    return eyedata


def simdir(data, TDir, TDur):
    """Simplify scanpaths based on angular relations between saccades (direction).

    Simplify consecutive saccades if the angle between them is smaller than the
    threshold TDir and the duration of the intermediate fixations is lower
    than threshold TDur.

    Parameters
    ---------
    data: array-like
        list of lists, output of gen_scanpath_structure
    TDir: float
        angle in degrees
    TDur: float
        time in seconds

    Returns
    ---------
    eyedata: list of lists, one iteration of direction based simplification
    """

    if len(data[3]) < 1:
        return data
    # if the scanpath is long enough
    else:
        i = 0
        j = 0
        # initialize empty lists
        sim_dur = []
        sim_x = []
        sim_y = []
        sim_lenx = []
        sim_leny = []
        sim_theta = []
        sim_len = []
        # while we don't run into index errors
        while i <= len(data[3]) - 1:
            if i < len(data[3]) - 1:
                # lets check angles
                v1 = [data[5][i], data[6][i]]
                v2 = [data[5][i + 1], data[6][i + 1]]
                angle = calcangle(v1, v2)
            else:
                # an angle of infinite size won't go into any further loop
                angle = float('inf')
            # if the angle is small and its not the last saccade
            if (angle < TDir) & (i < len(data[3]) - 1):
                # if the fixation duration is short:
                if data[2][i + 1] < TDur:
                    # if the fixation durations are short:
                    # calculate the sum of local vectors
                    v_x = data[5][i] + data[5][i + 1]
                    v_y = data[6][i] + data[6][i + 1]
                    rho, theta = cart2pol(v_x, v_y)
                    # save them in the new vectors
                    sim_lenx.insert(j, v_x)
                    sim_leny.insert(j, v_y)
                    sim_x.insert(j, data[3][i])
                    sim_y.insert(j, data[4][i])
                    sim_theta.insert(j, theta)
                    sim_len.insert(j, rho)
                    # add the fixation duration
                    sim_dur.insert(j, data[2][i])
                    i += 2
                    j += 1
                else:
                    # insert original data in new list
                    sim_lenx.insert(j, data[5][i])
                    sim_leny.insert(j, data[6][i])
                    sim_x.insert(j, data[3][i])
                    sim_y.insert(j, data[4][i])
                    sim_theta.insert(j, data[7][i])
                    sim_len.insert(j, data[8][i])
                    sim_dur.insert(j, data[2][i])
                    j += 1
                    i += 1
            # elif the angle is small, but its the last saccade:
            elif (angle < TDir) & (i == len(data[3]) - 1):
                # if the fixation duration is short:
                if data[2][i + 1] < TDur:
                    # calculate sum of local vectors
                    v_x = data[5][i - 2] + data[5][i - 1]
                    v_y = data[6][i - 2] + data[6][i - 1]
                    rho, theta = cart2pol(v_x, v_y)
                    # save them in new vectors
                    sim_lenx[j - 1] = v_x
                    sim_leny[j - 1] = v_y
                    sim_theta[j - 1] = theta
                    sim_len[j - 1] = rho
                    sim_dur.insert(j, data[2][-1] + (data[2][i] / 2))
                    j -= 1
                    i += 1
                # if fixation duration is long:
                else:
                    # insert original data in new list
                    sim_lenx.insert(j, data[5][i])
                    sim_leny.insert(j, data[6][i])
                    sim_x.insert(j, data[3][i])
                    sim_y.insert(j, data[4][i])
                    sim_theta.insert(j, data[7][i])
                    sim_len.insert(j, data[8][i])
                    sim_dur.insert(j, data[2][i])
                    i += 1
                    j += 1
            # else (the angle is too large
            else:
                # insert original data in new list
                sim_lenx.insert(j, data[5][i])
                sim_leny.insert(j, data[6][i])
                sim_x.insert(j, data[3][i])
                sim_y.insert(j, data[4][i])
                sim_theta.insert(j, data[7][i])
                sim_len.insert(j, data[8][i])
                sim_dur.insert(j, data[2][i])
                i += 1
                j += 1
    # now append the last fixation duration
    sim_dur.append(data[2][-1])
    eyedata = [
               [],
               [],
               sim_dur,
               sim_x,
               sim_y,
               sim_lenx,
               sim_leny,
               sim_theta,
               sim_len
               ]
    return eyedata


def simplify_scanpath(data, TAmp, TDir, TDur):
    """Simplify scanpaths until no further simplification is possible.

    Loops over simplification functions simdir and simlen until no
    further simplification of the scanpath is possible.

    Parameters
    -----------
    data: list of lists,
        output of gen_scanpath_structure
    TAmp: float
        length in px
    TDir: float
        angle in degrees
    TDur: float
        duration in seconds

    Returns
    -----------
    eyedata: list of lists, simplified vector-based scanpath representation
    """
    looptime = 0
    while True:
        data = simdir(data, TDir, TDur)
        data = simlen(data, TAmp, TDur)
        looptime += 1
        if looptime == len(data[2]):
            return data
            break


def cal_vectordifferences(data1, data2):
    """Create matrix of vector-length differences of all vector pairs

    Create M, a Matrix with all possible saccade-length differences between
    saccade pairs.

    Parameters
    -----------
    data1, data2: list of lists, vector-based scanpath representations

    Returns
    -----------
    M: array-like
        Matrix of vector length differences

    """
    # take length in x and y direction of both scanpaths
    x1 = np.asarray(data1[5])
    x2 = np.asarray(data2[5])
    y1 = np.asarray(data1[6])
    y2 = np.asarray(data2[6])
    # initialize empty lists M and row, will become matrix to store sacc-length
    # pairings
    M = []
    row = []
    # calculate saccade length differences, vectorized
    for i in range(0, len(x1)):
        x_diff = abs(x1[i] * np.ones(len(x2)) - x2)
        y_diff = abs(y1[i] * np.ones(len(y2)) - y2)
        # calc final length from x and y lengths, append, stack into matrix M
        row.append(np.asarray(np.sqrt(x_diff ** 2 + y_diff ** 2)))
        M = np.stack(row)
    return M


def createdirectedgraph(szM, M, M_assignment):
    """Create a directed graph
    The data structure of the result is a dicitionary within a dictionary
    such as
    weightedGraph = {0 : {1:259.55, 15:48.19, 16:351.95},
    1 : {2:249.354, 16:351.951, 17:108.97},
    2 : {3:553.30, 17:108.97, 18:341.78}, ...}

    Parameters
    -----------
    szM: list
        shape of matrix M
    M: array-like
        matrix of vector length differences
    M_assignment: array-like
        Matrix, arranged with values from 0 to number of entries in M

    Returns
    -----------
    weighted graph: dict
        Dictionary within a dictionary pairing weights (distances) with
        node-pairings

    """

    # initialize dictionary for neighbouring vertices and edge weights
    adjacent = {}
    weight = {}
    # loop through every node rowwise
    for i in range(0, szM[0]):
        # loop through every node columnwise
        for j in range(0, szM[1]):
            currentNode = i * szM[1] + j
            # if in the last (bottom) row, only go right
            if (i == szM[0] - 1) & (j < szM[1] - 1):
                adjacent[M_assignment[i, j]] = [currentNode + 1]
                weight[M_assignment[i, j]] = [M[i, j + 1]]
            # if in the last (rightmost) column, only go down
            elif (i < szM[0] - 1) & (j == szM[1] - 1):
                adjacent[M_assignment[i, j]] = [currentNode + szM[1]]
                weight[M_assignment[i, j]] = [M[i + 1, j]]
            # if in the last (bottom-right) vertex, do not move any further
            elif (i == szM[0] - 1) & (j == szM[1] - 1):
                adjacent[M_assignment[i, j]] = [currentNode]
                weight[M_assignment[i, j]] = [0]
            # anywhere else, move right, down and down-right.
            else:
                adjacent[M_assignment[i, j]] = [currentNode + 1,
                                                currentNode + szM[1],
                                                currentNode + szM[1] + 1]
                weight[M_assignment[i, j]] = [M[i, j + 1],
                                              M[i + 1, j],
                                              M[i + 1, j + 1]]
    # create list of all Nodes
    # Nodes = np.hstack(list(adjacent.values()))
    # create list of all associated distances (=weights)
    # Distances = np.hstack(list(weight.values()))
    # create ascending list ranging from first to last node
    Startnodes = range(0, szM[0] * szM[1])
    # initialize list with adjacent nodes and their weights
    weightedEdges = []
    # zip Nodes and weights
    for i in range(0, len(adjacent)):
        weightedEdges.append(list(zip(list(adjacent.values())[i],
                                      list(weight.values())[i])))
    # initialize final dictionary
    weightedGraph = {}
    # zip Startnodes together with Nodes-Weights, result is a nested dict
    for i in range(0, len(weightedEdges)):
        weightedGraph[Startnodes[i]] = dict(weightedEdges[i])
    return weightedGraph


def dijkstra(weightedGraph, start, end):
    """Implementation of Dijkstra algorithm
    Use the dijkstra algorithm to find the shortest path through a directed
    graph (weightedGraph) from start to end.

    Parameters
    ------------
    weightedGraph: dict
        dictionary within a dictionary pairing weights (distances) with
        node-pairings
    start: int
        starting point of path, should be 0
    end: int
        end point of path, should be (n, m) of Matrix M

    Returns
    ------------
    path: array-like
        array of indices of the shortest path, i.e. best-fitting saccade pairs
    dist: float
        sum of weights

    References
    ------------
    Dijkstra, E. W. (1959). A note on two problems in connexion with graphs.
    Numerische mathematik, 1(1), 269-271.
    """

    # initialize empty dictionary to hold distances
    dist = {}
    # inialize list of vertices in the path to current vertex (predecessors)
    pred = {}
    # where do I need to go?
    to_assess = weightedGraph.keys()
    for node in weightedGraph:
        # set inital distances to infinity
        dist[node] = float('inf')
        # no node has any predecessors yet
        pred[node] = None
    # initialize list to be filled with final distances(weights) of nodes
    sp_set = []
    # the starting node get a weight of 0 to make sure to start there
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


def cal_angulardifference(data1, data2, path, M_assignment):
    """Calculate angular similarity of two scanpaths.

    Parameters
    ------------
    data1: array-like
        list of lists, contains vector-based scanpath representation of the
        first scanpath
    data2: array-like
        list of lists, contains vector-based scanpath representation of the
        second scanpath
    path: array-like
        array of indices for the best-fitting saccade pairings between scan-
        paths
    M_assignment: array-like
         Matrix, arranged with values from 0 to number of entries in M, the
         matrix of vector length similarities

    Returns
    ------------
    anglediff: array of floats
        array of angular differences between pairs of saccades of two scanpaths

    """
    # get the angle between saccades from the scanpaths
    theta1 = data1[7]
    theta2 = data2[7]
    # initialize list to hold individual angle differences
    anglediff = []
    # calculate angular differences between the saccades along specified path
    for k in range(0, len(path)):
        # which saccade indices correspond to path?
        i, j = np.where(M_assignment == path[k])
        # extract the angle
        spT = [theta1[np.asscalar(i)], theta2[np.asscalar(j)]]
        for t in range(0, len(spT)):
            # get results in range -pi, pi
            if spT[t] < 0:
                spT[t] = math.pi + (math.pi + spT[t])
        spT = abs(spT[0] - spT[1])
        if spT > math.pi:
            spT = 2 * math.pi - spT
        anglediff.append(spT)
    return anglediff


def cal_durationdifference(data1, data2, path, M_assignment):
    """Calculate similarity of two scanpaths fixation durations.

    Parameters
    -----------

    :param data1: array-like
        list of lists, contains vector-based scanpath representation of the
        first scanpath
    data2: array-like
        list of lists, contains vector-based scanpath representation of the
        second scanpath
    :param path: array-like
        array of indices for the best-fitting saccade pairings between scan-
        paths
    :param M_assignment: array-like
         Matrix, arranged with values from 0 to number of entries in M, the
         matrix of vector length similarities

    Returns
    ----------
    durdiff: array of floats,
        array of fixation duration differences between pairs of saccades from
        two scanpaths

    """
    # get the duration of fixations in the scanpath
    dur1 = data1[2]
    dur2 = data2[2]
    # initialize list to hold individual duration differences
    durdiff = []
    # calculation fixation duration differences between saccades along path
    for k in range(0, len(path)):
        # which saccade indices correspond to path?
        i, j = np.where(M_assignment == path[k])
        maxlist = [dur1[np.asscalar(i)], dur2[np.asscalar(j)]]
        # compute abs. duration diff, normalize by largest duration in pair
        durdiff.append(abs(dur1[np.asscalar(i)] -
                           dur2[np.asscalar(j)]) / abs(max(maxlist)))
    return durdiff


def cal_lengthdifference(data1, data2, path, M_assignment):
    """Calculate length similarity of two scanpaths.

    Parameters
    ------------
    data1: array-like
        list of lists, contains vector-based scanpath representation of the
        first scanpath
    data2: array-like
        list of lists, contains vector-based scanpath representation of the
        second scanpath
    path: array-like
        array of indices for the best-fitting saccade pairings between scan-
        paths
    M_assignment: array-like
         Matrix, arranged with values from 0 to number of entries in M, the
         matrix of vector length similarities

    Returns
    ------------
    lendiff: array of floats
        array of length difference between pairs of saccades of two scanpaths

    """
    # get the saccade lengths rho
    len1 = np.asarray(data1[8])
    len2 = np.asarray(data2[8])
    # initialize list to hold individual length differences
    lendiff = []
    # calculate length differences between saccades along path
    for k in range(0, len(path)):
        i, j = np.where(M_assignment == path[k])
        lendiff.append(abs(len1[i] - len2[j]))
    return lendiff


def cal_positiondifference(data1, data2, path, M_assignment):
    """Calculate position similarity of two scanpaths.

    Parameters
    ------------
    data1: array-like
        list of lists, contains vector-based scanpath representation of the
        first scanpath
    data2: array-like
        list of lists, contains vector-based scanpath representation of the
        second scanpath
    path: array-like
        array of indices for the best-fitting saccade pairings between scan-
        paths
    M_assignment: array-like
         Matrix, arranged with values from 0 to number of entries in M, the
         matrix of vector length similarities

    Returns
    ------------
    posdiff: array of floats
        array of position differences between pairs of saccades
        of two scanpaths

    """
    # get the x and y coordinates of points between saccades
    x1 = np.asarray(data1[3])
    x2 = np.asarray(data2[3])
    y1 = np.asarray(data1[4])
    y2 = np.asarray(data2[4])
    # initialize list to hold individual position differences
    posdiff = []
    # calculate position differences along path
    for k in range(0, len(path)):
        i, j = np.where(M_assignment == path[k])
        posdiff.append(math.sqrt((x1[np.asscalar(i)] - x2[np.asscalar(j)]) ** 2 +
                                 (y1[np.asscalar(i)] - y2[np.asscalar(j)]) ** 2))
    return posdiff


def cal_vectordifferencealongpath(data1, data2, path, M_assignment):
    """Calculate vector similarity of two scanpaths.

    Parameters
    ------------
    data1: array-like
        list of lists, contains vector-based scanpath representation of the
        first scanpath
    data2: array-like
        list of lists, contains vector-based scanpath representation of the
        second scanpath
    path: array-like
        array of indices for the best-fitting saccade pairings between scan-
        paths
    M_assignment: array-like
         Matrix, arranged with values from 0 to number of entries in M, the
         matrix of vector length similarities

    Returns
    ------------
    vectordiff: array of floats
        array of vector differences between pairs of saccades of two scanpaths

    """
    # get the saccade lengths in x and y direction of both scanpaths
    x1 = np.asarray(data1[5])
    x2 = np.asarray(data2[5])
    y1 = np.asarray(data1[6])
    y2 = np.asarray(data2[6])
    # initialize list to hold individual vector differences
    vectordiff = []
    # calculate vector differences along path
    for k in range(0, len(path)):
        i, j = np.where(M_assignment == path[k])
        vectordiff.append(np.sqrt((x1[np.asscalar(i)] - x2[np.asscalar(j)]) ** 2 +
                                  (y1[np.asscalar(i)] - y2[np.asscalar(j)]) ** 2))
    return vectordiff


def getunnormalised(data1, data2, path, M_assignment):
    """Calculate unnormalised similarity measures.

    Calls the five functions to create unnormalised similarity measures for
    each of the five similarity dimensions. Takes the median of the resulting
    similarity values per array.

    Parameters
    ------------
    data1: array-like
        list of lists, contains vector-based scanpath representation of the
        first scanpath
    data2: array-like
        list of lists, contains vector-based scanpath representation of the
        second scanpath
    path: array-like
        array of indices for the best-fitting saccade pairings between scan-
        paths
    M_assignment: array-like
         Matrix, arranged with values from 0 to number of entries in M, the
         matrix of vector length similarities

    Returns
    -----------
    unnormalised: array
        array of unnormalised similarity measures on five dimensions

    Examples
    -----------
    >>> unorm_res = getunnormalised(scanpath_rep1, scanpath_rep2, path, M_assignment)
    """
    args = data1, data2, path, M_assignment
    VecSim = np.median(cal_vectordifferencealongpath(*args))
    DirSim = np.median(cal_angulardifference(*args))
    LenSim = np.median(cal_lengthdifference(*args))
    PosSim = np.median(cal_positiondifference(*args))
    DurSim = np.median(cal_durationdifference(*args))
    unnormalised = [VecSim, DirSim, LenSim, PosSim, DurSim]
    return unnormalised


def normaliseresults(unnormalised, sz=[1280, 720]):
    """Normalize similarity measures.

    Vector similarity is normalised against two times screen diagonal,
    the maximum
    theoretical distance.
    Direction similarity is normalised against pi.
    Length Similarity is normalised against screen diagonal.
    Position Similarity and Duration Similarity are already normalised.

    Parameters
    ------------
    unnormalised: array
        array of unnormalised similarity measures,
        output of getunnormalised()

    Returns
    ------------
    normalresults: array
        array of normalised similarity measures

    Examples
    ------------
    >>> normal_res = normaliseresults(unnormalised, sz = [1280, 720])
    """
    # normalize vector similarity against two times screen diagonal, the maximum
    # theoretical distance
    VectorSimilarity = 1 - unnormalised[0] / (2 * math.sqrt(sz[0] ** 2 + sz[1] ** 2))
    # normalize against pi
    DirectionSimilarity = 1 - unnormalised[1] / math.pi
    # normalize against screen diagonal
    LengthSimilarity = 1 - unnormalised[2] / math.sqrt(sz[0] ** 2 + sz[1] ** 2)
    PositionSimilarity = 1 - unnormalised[3] / math.sqrt(sz[0] ** 2 + sz[1] ** 2)
    # no normalisazion necessary, already done
    DurationSimilarity = 1 - unnormalised[4]
    normalresults = [VectorSimilarity, DirectionSimilarity, LengthSimilarity,
                     PositionSimilarity, DurationSimilarity]
    return normalresults


def docomparison(fixation_vectors1,
                 fixation_vectors2,
                 sz,
                 grouping,
                 TDir,
                 TDur,
                 TAmp
                 ):
    """Compare two scanpaths on five similarity dimensions.

     Parameters
     ------------

    fixation_vectors1: array-like
        n x 3 fixation vector of one scanpath
    fixation_vectors2: array-like
        n x 3 fixation vector of one scanpath
    sz: list
        screen dimensions in px.
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
        array of 5 scanpath similarity measures. Vector (Shape), Direction
    (Angle), Length, Position, and Duration. 1 means absolute similarity, 0 means
    lowest similarity possible.

    Examples
    ------------
    >>> results = docomparison(fix_1, fix_2, sz = [1280, 720], grouping = True, TDir = 45.0, TDur = 0.05, TAmp = 150)
    >>> print(results)
    >>> [[0.95075847681364678, 0.95637548674423822, 0.94082367355291008, 0.94491164030498609, 0.78260869565217384]]
    """
    # initialize result vector
    scanpathcomparisons = []
    # check if fixation vectors/scanpaths are long enough
    if (len(fixation_vectors1) >= 3) & (len(fixation_vectors2) >= 3):
        subj1 = gen_scanpath_structure(fixation_vectors1)
        subj2 = gen_scanpath_structure(fixation_vectors2)
        if grouping:
            subj1 = simplify_scanpath(subj1, TAmp, TDir, TDur)
            subj2 = simplify_scanpath(subj2, TAmp, TDir, TDur)
        M = cal_vectordifferences(subj1, subj2)
        szM = np.shape(M)
        M_assignment = np.arange(szM[0] * szM[1]).reshape(szM[0], szM[1])
        weightedGraph = createdirectedgraph(szM, M, M_assignment)
        path, dist = dijkstra(weightedGraph, 0, szM[0] * szM[1] - 1)
        unnormalised = getunnormalised(subj1, subj2, path, M_assignment)
        normal = normaliseresults(unnormalised, sz)
        scanpathcomparisons.append(normal)
    # return nan as result if at least one scanpath it too short
    else:
        scanpathcomparisons.append(np.repeat(np.nan, 5))
    return scanpathcomparisons

def main():
    print('Im executing main().')
    result = docomparison(data1,
                          data2,
                          sz,
                          grouping,
                          TDir,
                          TDur,
                          TAmp)
    print('Vector similarity = ', result[0][0])
    print('Direction similarity = ', result[0][1])
    print('Length similarity = ', result[0][2])
    print('Position similarity = ', result[0][3])
    print('Duration similarity = ', result[0][4])

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    # define arguments
    parser.add_argument('-i', '--input1', nargs='+', help='Input1: eyemovement data of the first subject',
                        metavar='PATH', required=True)
    parser.add_argument('-j', '--input2', nargs='+', help='Input2: eyemovement data of the second subject',
                        metavar='PATH', required=True)
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

    args = parser.parse_args()

    # read in data
    data1 = np.recfromcsv(args.input1[0], delimiter='\t',
                          dtype={'names': ('start_x', 'start_y', 'duration'), 'formats': ('f8', 'f8', 'f8')})
    data2 = np.recfromcsv(args.input2[0], delimiter='\t',
                          dtype={'names': ('start_x', 'start_y', 'duration'), 'formats': ('f8', 'f8', 'f8')})

    TDir = args.direction_threshold
    TAmp = args.amplitude_threshold
    TDur = args.duration_threshold
    sz = args.screensize

    if (TDir != 0) and (TAmp != 0):
        grouping = True
        print(
            'Scanpath comparison is done with grouping saccades shorter than {}px and with an angle smaller than {} degrees'
            ' if consecutive fixation are shorter than {} seconds.'.format(TAmp, TDir, TDur))
    else:
        grouping = False
        print('Scanpath comparison is done without any grouping')

    # execution
    main()
#    result = docomparison(data1,
#                          data2,
#                          sz,
#                          grouping,
#                          TDir,
#                          TDur,
#                          TAmp)

#    print('Vector similarity = ', result[0][0])
#    print('Direction similarity = ', result[0][1])
#    print('Length similarity = ', result[0][2])
#    print('Position similarity = ', result[0][3])
#    print('Duration similarity = ', result[0][4])
