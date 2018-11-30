#!/usr/bin/python3

#testdata in /home/data/psyinf/scratch/studyforrest
#tested until FixationChunks
#run in scratch/studyforrest-data-eyemovementlabels
#npdata = np.recfromcsv("sub-01/sub-01_task-movie_run-3_events.tsv",
#            delimiter='\t',
#            dtype={'names':('onset', 'duration', 'label', 'start_x', 'start_y',
#            'end_x', 'end_y', 'amp', 'peak_vel', 'med_vel', 'avg_vel'),
#            'formats':('f8', 'f8', 'U10', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',
#            'f8', 'f8')})

#npdata2 = np.recfromcsv("sub-02/sub-02_task-movie_run-3_events.tsv",
#            delimiter='\t',
#            dtype={'names':('onset', 'duration', 'label', 'start_x', 'start_y',
#            'end_x', 'end_y', 'amp', 'peak_vel', 'med_vel', 'avg_vel'),
#            'formats':('f8', 'f8', 'U10', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',
#            'f8', 'f8')})

#shots = pd.read_csv("../../../../../home/adina/MA_test/inputs/studyforrest-data-annotations/segments/avmovie/locations_run-3_events.tsv", sep = '\t')



import numpy as np
import pandas as pd
import math
from bisect import bisect_right
from bisect import bisect_left
import os
import sys, os
sys.path.insert(0, os.path.abspath('./'))
import MultiMatch.MultiMatch_pure as Mp

def cart2pol(x, y):
    '''transform cartesian into polar coordinates. Returns rho (length from 0,0)
    and theta (angle).'''
    rho = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y,x)
    return(rho, theta)

def calcangle(x1, x2):
    '''Calculate and return the angle between to vectors (saccades).'''
    angle = math.degrees(
            math.acos(
            np.dot(x1, x2)/(np.linalg.norm(x1)*np.linalg.norm(x2))))
    return angle

#Functions ported or adapted from MultiMatch
def generateStructureArrayScanpath(data):
    '''Take an n x 3 fixation vector (start_x, start_y, duration) in the form of 
    of a record array and transform it into appropriate vectorbased scanpath
    representation. Indices are as follows:
    0: fixation_x
    1: fixation_y
    2: fixation_dur
    3: saccade_x
    4: saccade_y
    5: saccade_lenx
    6: saccade_leny
    7: saccade_theta
    8: saccade_rho'''
    #initialize empty lists
    fixation_x = []
    fixation_y = []
    fixation_dur = []
    saccade_x = []
    saccade_y = []
    saccade_lenx = []
    saccade_leny = []
    saccade_theta = []
    saccade_rho = []
    #get the number of rows
    length=np.shape(data)[0]
    #keep coordinates and durations of fixations
    for i in range(0, length):
        fixation_x.append(data[i]['start_x'])
        fixation_y.append(data[i]['start_y'])
        fixation_dur.append(data[i]['duration'])
    #fixations as start coordinates for saccades (ignores PSOs and the like)
    for i in range(0, length-1):
        saccade_x.append(data[i]['start_x'])
        saccade_y.append(data[i]['start_y'])
    #calculate saccade length and angle
    for i in range(1, length):
        saccade_lenx.append(fixation_x[i] - saccade_x[i-1])
        saccade_leny.append(fixation_y[i] - saccade_y[i-1])
        rho, theta = cart2pol(saccade_lenx[i-1], saccade_leny[i-1])
        saccade_rho.append(rho)
        saccade_theta.append(theta)
    #append everything. Eyedata is a list of lists.
    eyedata = [fixation_x, fixation_y, fixation_dur, saccade_x, saccade_y,
              saccade_lenx, saccade_leny, saccade_theta, saccade_rho]
    return eyedata

def calVectordifferences(data1, data2):
    '''create M, a Matrix with all possible saccade-length differences between
    saccade pairs. Takes two scanpaths from generateSrtuctureArray() as input.'''
    #take length in x and y direction of both scanpaths
    x1 = np.asarray(data1[5])
    x2 = np.asarray(data2[5])
    y1 = np.asarray(data1[6])
    y2 = np.asarray(data2[6])
    #initialize empty lists M and row, will become matrix to store sacc-length
    #pairings
    M = []
    row = []
    #calculate saccade length differences, vectorized
    for i in range(0, len(x1)):
        x_diff = abs(x1[i] * np.ones(len(x2)) - x2)
        y_diff = abs(y1[i] * np.ones(len(y2)) - y2)
        #calc final length from x and y lengths, append, stack into matrix M
        row.append(np.asarray(np.sqrt(x_diff**2 + y_diff**2)))
        M = np.stack(row)
    return M


def createdirectedgraph(szM, M, M_assignment):
    '''create a directed graph. The data structure of the result is a
    dicitionary within a dictionary as in
    https://stackoverflow.com/questions/22897209/dijkstras-algorithm-in-python:
    weightedGraph = {0 : {1:259.55, 15:48.19, 16:351.95}, 1 : {2:249.354, 16:351.951,
    17:108.97}, 2 : {3:553.30, 17:108.97, 18:341.78}, ...}
    szM = shape of M
    M = Matrix with all saccade-length difference pairings
    M_assignment = Matrix, aranged with values from 0 to end if szM[0]*szM[1]
    entries'''
    #initialize dictionary for neighbouring vertices and edge weights
    adjacent = {}
    weight = {}
    #loop through every node rowwise
    for i in range(0, szM[0]):
        #loop through every node columnwise
        for j in range(0, szM[1]):
            currentNode = i * szM[1]+j
            #if in the last (bottom) row, only go right
            if (i == szM[0]-1) & (j < szM[1]-1):
                adjacent[M_assignment[i,j]] = [currentNode + 1]
                weight[M_assignment[i,j]] = [M[i, j+1]]
            #if in the last (rightmost) column, only go down
            elif (i < szM[0]-1) & (j ==szM[1]-1):
                adjacent[M_assignment[i,j]] = [currentNode + szM[1]]
                weight[M_assignment[i, j]] = [M[i + 1, j]]
            #if in the last (bottom-right) vertex, do not move any further
            elif (i == szM[0]-1) & (j == szM[1]-1):
                adjacent[M_assignment[i, j]] = [currentNode]
                weight[M_assignment[i, j]] = [0]
            #anywhere else, move right, down and down-right.
            else:
                adjacent[M_assignment[i, j]] = [currentNode + 1, currentNode + szM[1], currentNode + szM[1]+1]
                weight[M_assignment[i, j]] = [M[i, j+1],M[i + 1, j], M[i + 1, j + 1]]
    #create list of all Nodes
    Nodes = np.hstack(list(adjacent.values()))
    #create list of all associated distances (=weights)
    Distances = np.hstack(list(weight.values()))
    #create ascending list ranging from first to last node
    Startnodes = range(0, szM[0] * szM[1])
    #initialize list with adjacent nodes and their weights
    weightedEdges = []
    #zip Nodes and weights
    for i in range(0, len(adjacent)):
        weightedEdges.append(list(zip(list(adjacent.values())[i], list(weight.values())[i])))
    #initialize final dictionary
    weightedGraph = {}
    #zip Startnodes together with Nodes-Weights, result is a nested dict
    for i in range(0, len(weightedEdges)):
        weightedGraph[Startnodes[i]] = dict(weightedEdges[i])
    return weightedGraph


def dijkstra(weightedGraph, start, end):
    '''use the dijkstra algorithm to find the shortest path through a directed
    graph (weightedGraph) from start to end. weightedGraph is the output from
    createdirectedgraph() and has {node: {adjacent:weight}} structure.
    Returns the path and the final distance.'''
    #initialize empty dictionary to hold distances
    dist = {}
    #inialize list of vertices in the path to the current vertex (predecessors)
    pred = {}
    #where do I need to go?
    to_assess = weightedGraph.keys()
    for node in weightedGraph:
        #set inital distances to infinity
        dist[node] = float('inf')
        #no node has any predecessors yet
        pred[node] = None
    #initialize list to be filled with final distances(weights) of nodes
    sp_set = []
    #the starting node get a weight of 0 to make sure to start there
    dist[start] = 0
    #continue the algorithm as long as there are still unexplored nodes
    while len(sp_set) < len(to_assess):
        still_in = {node : dist[node] for node in [node for node in to_assess if
        node not in sp_set]}
        #find adjacent node with minimal weight and append to sp_set
        closest = min(still_in, key = dist.get)
        sp_set.append(closest)
        for node in weightedGraph[closest]:
            if dist[node] > dist[closest] + weightedGraph[closest][node]:
                dist[node] = dist[closest] + weightedGraph[closest][node]
                pred[node] = closest
    #append endnode to list path
    path = [end]
    #append contents of pred in reversed order to path
    while start not in path:
        path.append(pred[path[-1]])
    #return path in reverse order (begin to end) and final distance
    return path[::-1], dist[end]

def calAngularDifference(data1, data2, path, M_assignment):
    '''calculate the angular similarity of two scanpaths. Returns a list of
    angle-differences of aligned saccades in range -pi, pi.
    data1, data2: two scanpaths, output from generateStructureArrayScanpaths()
    path: shortest path (lowest distance) to align the two scanpaths
    M_assigment: Matrix of the size of all scanpath pairs from data1 and data2'''
    #get the angle between saccades from the scanpaths
    theta1 = data1[7]
    theta2 = data2[7]
    #initialize list to hold individual angle differences
    anglediff = []
    #calculate the angular differences between the saccades along specified path
    for k in range(0, len(path)):
        #which saccade indices correspond to path?
        i, j = np.where(M_assignment == path[k])
        #extract the angle
        spT = [theta1[np.asscalar(i)], theta2[np.asscalar(j)]]
        for t in range(0, len(spT)):
            #get results in range -pi, pi
            if spT[t] < 0:
                spT[t] = math.pi + (math.pi + spT[t])
        spT = abs(spT[0]-spT[1])
        if spT > math.pi:
            spT= 2 * math.pi - spT
        anglediff.append(spT)
    return anglediff

def calDurationDifference(data1, data2, path, M_assignment):
    '''calculate the duration similarity of two scanpaths. Returns a list of
    absolute duration-differences of aligned saccades, scaled by largest
    duration in the compared pair. Should be non-negative.
    data1, data2: two scanpaths, output from generateStructureArrayScanpaths()
    path: shortest path (lowest distance) to align the two scanpaths
    M_assigment: Matrix of the size of all scanpath pairs from data1 and data2'''
    #get the duration of fixations in the scanpath
    dur1 = data1[2]
    dur2 = data2[2]
    #initialize list to hold individual duration differences
    durdiff = []
    #calculation fixation duration differences between saccades along path
    for k in range(0, len(path)):
        #which saccade indices correspond to path?
        i, j = np.where(M_assignment == path[k])
        maxlist = [dur1[np.asscalar(i)], dur2[np.asscalar(j)]]
        #compute abs. duration difference, normalize by largest duration in pair
        durdiff.append(abs(dur1[np.asscalar(i)] -
            dur2[np.asscalar(j)])/abs(max(maxlist)))
    return durdiff

def calLengthDifference(data1, data2, path, M_assignment):
    '''calculate the length similarity of two scanpaths. Returns a list of
    absolute length differences of aligned saccades.
    data1, data2: two scanpaths, output from generateStructureArrayScanpaths()
    path: shortest path (lowest distance) to align the two scanpaths
    M_assigment: Matrix of the size of all scanpath pairs from data1 and data2'''
    #get the saccade lengths rho
    len1 = np.asarray(data1[8])
    len2 = np.asarray(data2[8])
    #initialize list to hold individual length differences
    lendiff = []
    #calculate length differences between saccades along path
    for k in range(0, len(path)):
        i, j = np.where(M_assignment == path[k])
        lendiff.append(abs(len1[i]-len2[j]))
    return lendiff
#TODO: this is a list of arrays.. change? [array([[118]), array([63.66]), ...]

def calPositionDifference(data1, data2, path, M_assignment):
    '''calculate the fixation position similarity of two scanpaths. Returns a
    list of absolute position differences of aligned fixations.
    data1, data2: two scanpaths, output from generateStructureArrayScanpaths()
    path: shortest path (lowest distance) to align the two scanpaths
    M_assigment: Matrix of the size of all scanpath pairs from data1 and data2'''
    #get the x and y coordinates of points between saccades
    x1 = np.asarray(data1[3])
    x2 = np.asarray(data2[3])
    y1 = np.asarray(data1[4])
    y2 = np.asarray(data2[4])
    #initialize list to hold individual position differences
    posdiff = []
    #calculate position differences along path
    for k in range(0, len(path)):
        i, j = np.where(M_assignment == path[k])
        posdiff.append(math.sqrt((x1[np.asscalar(i)] - x2[np.asscalar(j)])**2 +
            (y1[np.asscalar(i)] - y2[np.asscalar(j)])**2))
    return posdiff

def calVectorDifferenceAlongPath(data1, data2, path, M_assignment):
    '''calculate the vector similarity of two scanpaths. Returns a
    list of absolute vector differences of aligned scanpaths.
    data1, data2: two scanpaths, output from generateStructureArrayScanpaths()
    path: shortest path (lowest distance) to align the two scanpaths
    M_assigment: Matrix of the size of all scanpath pairs from data1 and data2'''
    #get the saccade lengths in x and y direction of both scanpaths
    x1 = np.asarray(data1[5])
    x2 = np.asarray(data2[5])
    y1 = np.asarray(data1[6])
    y2 = np.asarray(data2[6])
    #initialize list to hold individual vector differences
    vectordiff = []
    #calculate vector differences along path
    for k in range(0, len(path)):
        i, j = np.where(M_assignment == path[k])
        vectordiff.append(np.sqrt((x1[np.asscalar(i)] - x2[np.asscalar(j)])**2 +
        (y1[np.asscalar(i)] - y2[np.asscalar(j)])**2))
    return vectordiff

def getunnormalised(data1, data2, path, M_assignment):
    '''calculate the unnormalised similarity measures for the five similarity
    dimensions. Return the median value of the resulting lists per dimension.
    data1, data2: two scanpaths, output from generateStructureArrayScanpaths()
    path: shortest path (lowest distance) to align the two scanpaths
    M_assignemtn: Matrix of the size of all scanpath pairs from data1 and
    data2'''
    VecSim = np.median(calVectorDifferenceAlongPath(data1, data2, path,
                M_assignment))
    DirSim = np.median(calAngularDifference(data1, data2, path, M_assignment))
    LenSim = np.median(calLengthDifference(data1, data2, path, M_assignment))
    PosSim = np.median(calPositionDifference(data1, data2, path, M_assignment))
    DurSim = np.median(calDurationDifference(data1, data2, path, M_assignment))
    unnormalised = [VecSim, DirSim, LenSim, PosSim, DurSim]
    return unnormalised

def normaliseresults(unnormalised, sz = [1280, 720]):
    #normalize vector similarity against two times screen diagonal, the maximum
    #theoretical distance
    VectorSimilarity = 1 - unnormalised[0] / (2 * math.sqrt(sz[0]**2 + sz[1]**2))
    #normalize against pi 
    DirectionSimilarity = 1 - unnormalised[1] / math.pi
    #normalize against screen diagonal
    LengthSimilarity = 1 - unnormalised[2] / math.sqrt(sz[0]**2 + sz[1]**2)
    PositionSimilarity = 1 - unnormalised[3] / math.sqrt(sz[0]**2 + sz[1]**2)
    #no normalisazion necessary, already done
    DurationSimilarity = 1 - unnormalised[4]
    normalresults = [VectorSimilarity, DirectionSimilarity, LengthSimilarity,
                    PositionSimilarity, DurationSimilarity]
    return normalresults


#Functions specifically for the data at hand

def takeClosestright(myList, myNumber):
    '''return the integer closest to 'myNumber' in an ordered list.(shamelessly
    stolen from
    https://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-
    closest-to-a-given-value)
    '''
    pos = bisect_right(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    after = myList[pos]
    return after

def takeClosestleft(myList, myNumber):
    '''return the integer closest to 'myNumber' in an ordered list that is lower
    than 'myNumber'''
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    return before

def createOnsets(data, dur):
    '''create onset times of all shots of 'dur' seconds of length
    data = dataframe, should be location annotation
    dur = duration in seconds'''
    onsets = []
    for index, row in data.iterrows():
        if row['duration'] >= dur:
            onsets.append(row['onset'])
    return onsets


def createOffsets(data, dur):
    '''create offset times of all shots of 'dur' seconds of length
    data = dataframe, should be location annotation
    dur = duration in seconds'''
    offsets = []
    for index, row in data.iterrows():
        if row['duration'] >= dur:
            #calculate end of shot by adding onset + duration, subtract an
            #epsilon to be really sure not to get into a cut
            offsets.append(row['onset']+row['duration']-0.03)
    return offsets


def createChunks(onsets, fixations, dur):
    '''Create and return start and end indices to chunk eye movement data into
    segments to compare scanpaths across.
    onsets = output from CreateOnsets()
    fixations = output from preprocess(). n x 4 np record array
    dur = durations of segments in seconds'''
    #initialize empty lists
    startidx, endidx = [], []
    for shotonset in onsets:
        start = takeClosestright(fixations['onset'], shotonset)
        startidx.append(np.where(fixations['onset']==start)[0].tolist())
        end = takeClosestright(fixations['onset'], shotonset + dur)
        endidx.append(np.where(fixations['onset']==end)[0].tolist())
    #flatten the nested lists
    startidx = [element for sublist in startidx for element in sublist]
    endidx = [element for sublist in endidx for element in sublist]
    return startidx, endidx


#createOffsetChunks changes to **last** seconds of shot!!

def createOffsetChunks(offsets, fixations, dur):
    '''Create and return start and end indices to chunk eye movement data into
    segments to compare scanpaths across.
    onsets = output from CreateOff(!!)sets()
    fixations = output from preprocess(). n x 4 np record array
    dur = durations of segments in seconds'''
    startidx, endidx = [], []
    for shotoffset in offsets:
        start = takeClosestright(fixations['onset'], shotoffset - dur)
        startidx.append(np.where(fixations['onset']==start)[0].tolist())
        end = takeClosestleft(fixations['onset'], shotoffset)
        endidx.append(np.where(fixations['onset'] == end)[0].tolist())
    #flatten the nested lists
    startidx = [element for sublist in startidx for element in sublist]
    endidx = [element for sublist in endidx for element in sublist]
    return startidx, endidx


def FixationsChunks(fixations, startid, endid):
    '''Chunk eye movement data into segments of approximate length to compute
    scanpath similarities. Output is returned as a n x 3 fixation vector.
    startid, endid = output from createChunks
    fixations = output from preprocess'''
    fixation_vector = []
    #slice fixation data according to indices, take columns start_x, start_y and
    #duration
    for idx in range(0, len(startid)):
        ind = fixations[startid[idx]:endid[idx]][["start_x", "start_y", "duration"]]
        fixation_vector.append(ind)
    return fixation_vector

def pursuits_to_fixations(npdata):
    '''this function takes a numpy record array from Asims eye-event detection
    algorithm. Start and end points of pursuits are transformed into a fixation,
    the pursuit movement can then be simplified as a saccade. The function
    returns a recordarray'''
    #initialize empty rec array of the same shape
    newdata=np.recarray((0,), dtype=[('onset', '<f8'),
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
    #reassemble rec array. split pursuits to use end and start as fixations later
    for i in range(0, len(npdata)):
        if npdata[i]['label']=='PURS':
            row_1 = npdata[i]
            row_1['duration'] = npdata[i]['duration']/2
            row_2 = row_1.copy()
            row_2['onset']+= row_2['duration']
            row_2['start_x'] = row_2['end_x']
            row_2['start_y'] = row_2['end_y']
            newdata = np.append(newdata, row_1)
            newdata = np.append(newdata, row_2)
        else:
            newdata = np.append(newdata, npdata[i])
    return newdata


def preprocess(data, sz=[1280, 720]):
    '''
    data = n x 11 rec array
    sz = screen measurements
    preprocesses a recordarray with eye events (from pursuits_to_fixations()
    function. Assumes the datafile is sorted by time. Will filter to include
    only Fixations and Pursuits-start/end-points, will check for out-of-bound
    gazes.
    Returns clean Fixation data with onset, start_x, start_y and duration '''
    #only fixations and pursuits
    Filterevents = data[np.logical_or(data['label']=='FIXA',
    data['label']=='PURS')]
    #within x coordinates?
    Filterxbounds = Filterevents[np.logical_and(Filterevents['start_x'] >= 0,
    Filterevents['start_x'] <= sz[0])]
    #within y coordinates?
    Filterybounds = Filterxbounds[np.logical_and(Filterxbounds['start_y'] >= 0,
    Filterxbounds['end_y'] <= sz[1])]
    #give me onset times, start_x, start_y and duration
    fixations = Filterybounds[["onset", "start_x", "start_y",
    "duration"]]
    return fixations
#TODO: this still gives warning: FutureWarning: Numpy has detected that you may
#be viewing or writing to an array returned by selecting multiple fields in a
#structured array.
#This code may break in numpy 1.13 because this will return a view instead of a
#copy -- see release notes for details



def longshot(shots, dur = 4.92):
    '''group movie shots without a cut together to obtain longer movie
    segments. This way, fewer but longer scanpaths are obtained. Example: use
    median shotlength of 4.92s.
    shots = dataframe, contains movie location annotation
    dur = length in seconds for movie shot
    '''
    #turn pandas dataframe shots into record array
    structshots = shots.to_records()
    i = 0
    while i < len(structshots):
        #break before running into index error
        if structshots[i] == structshots[-1]:
            break
        else:
            if (structshots[i]['duration'] < dur) & \
            (structshots[i+1]['duration'] < dur) & \
            (structshots[i]['locale'] == structshots[i+1]['locale']):
                #add durations together and delete second row
                structshots[i]['duration'] += structshots[i+1]['duration']
                structshots = np.delete(structshots, i+1, 0)
            else:
                i += 1
    aggregated = pd.DataFrame({'onset':structshots['onset'].tolist(),
    'duration': structshots['duration'].tolist()}, columns = ['onset',
    'duration'])
    return aggregated


def doComparison(shots, data1, data2, sz = [1280, 720], dur = 5):
    '''
    Compare the scanpaths evoked by the same moviesegment of two subjects.
    Return a vector of five similarity measures: Vector (Shape), Direction
    (Angle), Length, Position, and Duration. 1 means absolute similarity, 0 mean
    lowest similarity possible.
    shots: pandas dataframe with locaction annotation
    data1: eyemovement data of subject 1
    data2: eyemovement data of subject 2
    sz: screen measurements. Default is 1280 x 720 px
    dur: shot duration over which scanpaths should be compared in seconds,
    default is 3 seconds.
    '''
    #initialize result vector
    scanpathcomparisons = []
    #transform pursuits into fixations
    newdata1 = pursuits_to_fixations(data1)
    newdata2 = pursuits_to_fixations(data2)
    #preprocess input files
    fixations1 = preprocess(newdata1, sz)
    fixations2 = preprocess(newdata2, sz)
    shots = longshot(shots, dur) #for longer shots, change dur

    #for the first seconds in a shot for comparison take the following code
    ###onset = createOnsets(shots, dur)
    ###startid1, endid1 = createChunks(onset, fixations1, dur)
    ###startid2, endid2 = createChunks(onset, fixations2, dur)

    #for the last seconds in a shot for comparison take this code instead
    onset = createOffsets(shots, dur)
    startid1, endid1 = createOffsetChunks(onset, fixations1, dur)
    startid2, endid2 = createOffsetChunks(onset, fixations2, dur)
    fixation_vectors1 = FixationsChunks(fixations1, startid1, endid1)
    fixation_vectors2 = FixationsChunks(fixations2, startid2, endid2)
    #save onset and duration times, if valid ones can be calculated
    onset_times = []
    exact_durations = []
    for i in range(0, len(startid1)):
        onset_time = fixations1[startid1[i]]['onset']
        onset_times.append(onset_time)
        exact_duration = fixations1[endid1[i]]['onset']-fixations1[startid1[i]]['onset']
    #capture negative durations for invalid scanpaths
        if exact_duration > 0:
            exact_durations.append(exact_duration)
        else:
            exact_durations.append(np.nan)
#    onset_times = fixations1[startid1]['onset']
#    exact_durations = fixations1[endid1]['onset']-fixations1[startid1]['onset']
    #loop over all fixation vectors/scanpaths and calculate similarity
    for i in range(0, len(onset)):
        #check if fixation vectors/scanpaths are long enough
        if (len(fixation_vectors1[i]) >= 3) & (len(fixation_vectors2[i]) >=3):
            subj1 = generateStructureArrayScanpath(fixation_vectors1[i])
            subj2 = generateStructureArrayScanpath(fixation_vectors2[i])
            M = calVectordifferences(subj1, subj2)
            szM = np.shape(M)
            M_assignment = np.arange(szM[0]*szM[1]).reshape(szM[0], szM[1])
            weightedGraph = createdirectedgraph(szM, M, M_assignment)
            path, dist = dijkstra(weightedGraph, 0, szM[0]*szM[1]-1)
            unnormalised = getunnormalised(subj1, subj2, path, M_assignment)
            normal = normaliseresults(unnormalised, sz)
            scanpathcomparisons.append(normal)
        #return nan as result if at least one scanpath it too short
        else:
            scanpathcomparisons.append(np.repeat(np.nan, 5))
    return scanpathcomparisons, onset_times, exact_durations


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    #define arguments
    parser.add_argument('-i', '--input1', nargs = '+', help = 'Input1: eyemovement data of the first subject', metavar = 'PATH', required = True)
    parser.add_argument('-j', '--input2', nargs = '+', help = 'Input2: eyemovement data of the second subject', metavar = 'PATH', required = True)
    parser.add_argument('-k', '--input3', help = 'Input3: location annotation of the movie segment', metavar = 'PATH', required = True)
    parser.add_argument('-o', '--output', help = 'Output: Specify path where output should be saved', metavar = 'PATH', required = True)
    parser.add_argument('-d', '--duration', help = 'duration: Specify the approx. time of a scanpath in seconds, i.e. 3.0. Note: Scanpaths are extracted within a shot, not across shots! Long durations will lead to only few scanpaths', type = float, default = 4.92)
    parser.add_argument('-ld', '--lduration', help = 'duration: group short shots in the same locale (i.e. no change of scenes between them) together for longer scanpaths', type=float, default=None)
    parser.add_argument('-di', '--direction_threshold', help='direction_threshold: for direction based grouping. If 0: no grouping will be performed', type = float, default=0.0)
    parser.add_argument('-am', '--amplitude_threshold', help='amplitude_threshold: for amplitude based grouping. If 0: no grouping will be performed', type = float, default=0.0)
    parser.add_argument('-du', '--duration_threshold', help='duration_threshold: for direction based grouping.', type = float, default=0.0)
    parser.add_argument('-sz', '--screensize', help='screensize: Resolution of screen in px, default is [720, 1280]', default = [720, 1280])
    parser.add_argument('-pos', '--position_offset', help='position_offset: if True, scanpaths of dur length stop at shotoffset (instead of starting at shotonset', default=False)

    args = parser.parse_args()

#Data read in

    shots = pd.read_csv(args.input3, sep = '\t')
    data1 = np.recfromcsv(args.input1[0],
            delimiter='\t',
            dtype={'names':('onset', 'duration', 'label', 'start_x', 'start_y',
            'end_x', 'end_y', 'amp', 'peak_vel', 'med_vel', 'avg_vel'),
            'formats':('f8', 'f8', 'U10', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',
            'f8', 'f8')})

    data2 = np.recfromcsv(args.input2[0],
            delimiter='\t',
            dtype={'names':('onset', 'duration', 'label', 'start_x', 'start_y',
            'end_x', 'end_y', 'amp', 'peak_vel', 'med_vel', 'avg_vel'),
            'formats':('f8', 'f8', 'U10', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',
            'f8', 'f8')})
    dur = args.duration
    TDir = args.direction_threshold
    TAmp = args.amplitude_threshold
    TDur = args.duration_threshold
    sz = args.screensize
    ldur = args.lduration
    offset = args.position_offset
    #derive simple boolean variable to
    if (TDir != 0) and (TAmp != 0):
        grouping = True
        print('Scanpath comparison is done with grouping saccades shorter than {}px and with an angle smaller than {}°'
              ' if consecutive fixation are shorter than {} seconds.'.format(TAmp, TDir, TDur))
    else:
        grouping = False
        print('Scanpath comparison is done without any grouping')

#Execution
#TODO: optional arguments need to go here as well (sz and dur)
#TODO: output format needs changing. Should be BIDS format. "Onset", "Duration"
#(obsolete?), "segment"
    segment, onset, duration = doComparison(shots, data1, data2)
#transform list into numpy array
    segmentfinal = np.array(segment)
#stack stuff together
    results = np.column_stack((onset, duration, segmentfinal))
#to record array
#    final = np.core.records.fromarrays(results.transpose(), names = "onset,
#    duration, vector_sim, direction_sim, length_sim, position_sim,
#    duration_sim", formats = 'f8, f8, f8, f8, f8, f8, f8')
#save
    if not os.path.isdir(os.path.dirname(args.output)):
        os.makedirs(os.path.dirname(args.output))
    np.savetxt(args.output,
                results,
                fmt = '%f\t%f\t%f\t%f\t%f\t%f\t%f',
                delimiter='\t',
                header=\
                "onset\tduration\tvector_sim\tdirection_sim\tlength_sim\tposition_sim\tduration_sim",
                comments = ''
                )



