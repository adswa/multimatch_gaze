#!/usr/bin/python3

"""
short script to transform remodnav
(https://github.com/psychoinformatics-de/remodnav) outputs from naturalistic
stimulation with hollywood movie forrest gump into fixation vectors
"""

import numpy as np
import pandas as pd
import math
from bisect import bisect_right
from bisect import bisect_left
import os.path as op


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


def createOnsets(data, dur):
    '''create onset times of all shots of 'dur' seconds of length
    data = dataframe, should be location annotation
    dur = duration in seconds'''
    onsets = []
    for index, row in data.iterrows():
        if row['duration'] >= dur:
            onsets.append(row['onset'])
    return onsets



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

def longshot(shots, dur):
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



def savefile(fixation_vector, output, header):
    newheader = ''.join([w+'\t' for w in header]).strip()
    np.savetxt(output, fixation_vector, delimiter='\t', comments='', header =  newheader)


def run(data1, shots, sz, dur, subname):
    newdata1 = pursuits_to_fixations(data1)
    fixations1 = preprocess(newdata1, sz)
    shots = longshot(shots, dur)
    onset = createOnsets(shots, dur)
    startid1, endid1 = createChunks(onset, fixations1, dur)
    fixation_vectors1 = FixationsChunks(fixations1, startid1, endid1)
    header = fixation_vectors1[0].dtype.names
    for i in range(0, len(onset)):
        output = args.output + '/fixvectors/segment_' + str(i) + '_' + subname + '.tsv'
        print('saving file', i, 'into', output)
        savefile(fixation_vectors1[i], output, header)



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    #define arguments
    parser.add_argument('-i', '--input1', nargs = '+', help = 'Input: eyemovement data of one subject', metavar = 'PATH', required = True)
    parser.add_argument('-o', '--output', help = 'Output: Specify path where output should be saved', metavar = 'PATH', required = True)
    parser.add_argument('-sz', '--screensize', help = 'Screensize: what are the dimensions of the screen the stimulus was displayed on in px, e.g. [1280, 720]', default=[1280, 720])
    parser.add_argument('-s', '--shots', help = 'Input3: location annotation of the movie segment', metavar = 'PATH', required = True)
    parser.add_argument('-d', '--duration', help = 'approximate duration of video segments to derive fixation vectors from', default=5.0)
    args = parser.parse_args()

#Data read in

    data1 = np.recfromcsv(args.input1[0],
            delimiter='\t',
            dtype={'names':('onset', 'duration', 'label', 'start_x', 'start_y',
            'end_x', 'end_y', 'amp', 'peak_vel', 'med_vel', 'avg_vel'),
            'formats':('f8', 'f8', 'U10', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',
            'f8', 'f8')})
    shots = pd.read_csv(args.shots, sep = '\t')

    subname = op.basename(args.input1[0]).split('_')[0]

    if args.screensize:
        sz = args.screensize
    else:
        sz = [1280, 720]

    if args.duration:
        dur = float(args.duration)
    else:
        dur = 5.0

    #run everything
    run(data1, shots, sz, dur, subname)


