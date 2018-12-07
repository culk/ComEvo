import networkx as nx
import snap
import os

import datetime
import time

import numpy as np
import snap

import pdb

import config
from graph import *
import itertools
from networkx import edge_betweenness_centrality

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

def generateNetworkxGraph(graph):
    networkxGraph = nx.MultiGraph()
    for edgeI in graph.Edges():
        edgeId = edgeI.GetId()
        srcId = edgeI.GetSrcNId()
        dstId = edgeI.GetDstNId()

        if not networkxGraph.has_node(srcId):
            networkxGraph.add_node(srcId)

        if not networkxGraph.has_node(dstId):
            networkxGraph.add_node(dstId)

        #Add edges and weights
        #for now only keep 1 edge between nodes even if it has multi-edges
        if not networkxGraph.has_edge(srcId, dstId):
            edgeweight = 1.0
            networkxGraph.add_edge(srcId, dstId, weight=edgeweight)

    return networkxGraph

def GirvanNewmanCommunityDetection(graph):
    initialComponents = nx.number_connected_components(graph)
    ncomp = initialComponents
    pdb.set_trace()
    while ncomp <= initialComponents:
        bw = nx.edge_betweenness_centrality(graph, weight='weight')    #edge betweenness for G
        #find the edge with max centrality
        print str(bw)

        max_ = max(bw.values())

        #find the edge with the highest centrality and remove all of them if there is more than one!
        for k, v in bw.iteritems():
            if float(v) == max_:
                graph.remove_edge(k[0],k[1])    #remove the central edge
        ncomp = nx.number_connected_components(graph)    #recalculate the no of components

def UpdateDeg(A, nodes):
    pdb.set_trace()
    deg_dict = {}
    n = len(nodes)  #len(A) ---> some ppl get issues when trying len() on sparse matrixes!
    B = A.sum(axis = 1)
    for i in range(n):
        if i in nodes:
            #deg_dict[nodes[i]] = B[i, 0]
            deg_dict[i] = B[i, 0]
    return deg_dict

#compute the modularity of current split
def GirvanNewmanGetModularity(G, deg_, m_):
    pdb.set_trace()
    New_A = nx.adjacency_matrix(G)
    New_deg = {}
    New_deg = UpdateDeg(New_A, G.nodes())
    #Let's compute the Q
    comps = nx.connected_components(G)    #list of components    
    print 'No of communities in decomposed G: %d' % nx.number_connected_components(G)
    Mod = 0    #Modularity of a given partitionning
    for c in comps:
        EWC = 0    #no of edges within a community
        RE = 0    #no of random edges
        for u in c:
            EWC += New_deg[u]
            RE += deg_[u]        #count the probability of a random edge
        Mod += ( float(EWC) - float(RE*RE)/float(2*m_) )
    Mod = Mod/float(2*m_)
    if _DEBUG_:
        print "Modularity: %f" % Mod
    return Mod

def most_central_edge(G):
    centrality = edge_betweenness_centrality(G, weight='weight')
    return max(centrality, key=centrality.get)

def main():
    graph = Graph('partial')
    graph.print_summary()

    graph_snapshots = get_graph_snapshots(graph.graph, 86400*30)
    print str(graph_snapshots)

    networkxGraph = generateNetworkxGraph(graph.graph)

    #G = path_graph(10)
    #comp = nx.algorithms.community.centrality.girvan_newman(networkxGraph)
    comp = nx.algorithms.community.centrality.girvan_newman(networkxGraph, most_valuable_edge=most_central_edge)
    #for communities in itertools.islice(comp, 3):
    #    print(tuple(sorted(c) for c in communities))

    communitues = tuple(sorted(c) for c in next(comp))
    print str(communitues)

    #print len(networkxGraph.edges())

    #print graph.graph.GetEdges()
    pdb.set_trace()
    n = networkxGraph.number_of_nodes()    #|V|
    A = nx.adjacency_matrix(networkxGraph)    #adjacenct matrix

    m_ = 0.0    #the weighted version for number of edges
    for i in range(0,n):
        for j in range(0,n):
            m_ += A[i,j]
    m_ = m_/2.0
    print "m: %f" % m_

    Orig_deg = {}
    Orig_deg = UpdateDeg(A, networkxGraph.nodes())

    GirvanNewmanCommunityDetection(networkxGraph)
    mod = GirvanNewmanGetModularity(networkxGraph, Orig_deg, m_)
    print "Modularity = " + str(mod)
main()