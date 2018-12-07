import datetime
import time
import os

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import snap
import sys

import config
import utilities

from graph import Graph
from utilities import *
import pdb


def get_graph_index(edgeTime, currentTime, timeQuantum):
    #integer division on purpose
    return (currentTime - edgeTime) / timeQuantum

def get_snapshot_timerange(snapShotIndex, currentTime, timeQuantum):
    #returns the epoch time range for the graph. Will be a tuple of start time and end time
    startTime = currentTime - (snapShotIndex*(timeQuantum + 1))
    endTime = currentTime - (snapShotIndex*(timeQuantum))
    return (startTime, endTime)

def get_graph_snapshots(graph, timeQuantum):
    """
    Arg - takes in the graph, timeQuantum should be in seconds
    Returns the list of graph snapshots with their epoch time ranges
    """
    startTime = 1217567877
    currentEpochTime = int(time.time())
    
    numberOfSnapshots = get_graph_index(startTime, currentEpochTime, timeQuantum) + 1 #Handle OBOB

    graphSnapshots = []
    graphSnapshotRanges = []

    for i in range(numberOfSnapshots):
        #Add empty multi-graphs
        newGraph = snap.TNEANet.New()
        newGraph.AddIntAttrE('time')
        graphSnapshots.append(newGraph)

    for edgeI in graph.Edges():
        edgeId = edgeI.GetId()
        srcId = edgeI.GetSrcNId()
        dstId = edgeI.GetDstNId()
        try:
            edgeTimeStamp = graph.GetIntAttrDatE(edgeI.GetId(), 'time')
            currentEdgeSnapshotIndex = get_graph_index(edgeTimeStamp, currentEpochTime, timeQuantum)
            graphSnapshot = graphSnapshots[currentEdgeSnapshotIndex]
            
            #Add the nodes
            if not graphSnapshot.IsNode(srcId):
                graphSnapshot.AddNode(srcId)

            if not graphSnapshot.IsNode(dstId):
                graphSnapshot.AddNode(dstId)

            #Add the edge and assign the timestamp as an attribute
            snapshotEdgeId = graphSnapshot.AddEdge(srcId, dstId)

            graphSnapshot.AddIntAttrDatE(snapshotEdgeId, edgeTimeStamp, 'time')

        except:
            #Edge Doesn't have timestamp!!! Handle Later
            raise

    nonEmptySnapshotsWithRanges = []
    #return only graphs that have >0 edges
    for i in range(len(graphSnapshots)): 
        g = graphSnapshots[i]
        if g.GetEdges() > 0:
            nonEmptySnapshotsWithRanges.append((get_snapshot_timerange(i, currentEpochTime, timeQuantum), g))

    """
    Debug Block - 

    totalEdges = 0
    for g in graphSnapshots:
        if g.GetEdges() > 0:
            for e in g.Edges():
                print "e is " + str(e.GetSrcNId()) + " -> " + str(e.GetDstNId()) + ", time = " + str(g.GetIntAttrDatE(e.GetId(), 'time')) + ", edges = " + str(g.GetEdges())
            raw_input()
            totalEdges += g.GetEdges()

    print "Total Edges = " + str(totalEdges)

    """

    """
    Debug block

    g = graphSnapshots[3744]

    for e in g.Edges():
        print "e is " + str(e.GetSrcNId()) + " -> " + str(e.GetDstNId()) + ", time = " + str(g.GetIntAttrDatE(e.GetId(), 'time')) + ", edges = " + str(g.GetEdges())

    """
    return nonEmptySnapshotsWithRanges

    
def main():
    graph = Graph('30days')
    graph.print_summary()
    graph.preprocess()
    graph.set_time_delta(86400*3)
    graph.calc_communities("leiden-algorithm", weight_fn=linear_fn, weighted=True)
    #graph.calc_communities("fastgreedy", weight_fn=None, weighted=False)
    #graph.get_conductance()

if __name__ == '__main__':
    main()
