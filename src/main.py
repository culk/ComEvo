import os
import snap
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import config
import pdb
import time
import datetime

def get_graph_index(edgeTime, currentTime, timeQuantum):
	#integer division on purpose
	return (currentTime - edgeTime) / timeQuantum

def get_graph_snapshots(graph, timeQuantum):
	"""
	Arg - takes in the graph, timeQuantum should be in seconds
	Returns the list of graph snapshots
	"""
	#
	startTime = 1217567877
	currentEpochTime = int(time.time())
	print currentEpochTime
	
	numberOfSnapshots = get_graph_index(1217567877, currentEpochTime, timeQuantum) + 1 #Handle OBOB

	graphSnapshots = []

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
			#edgeTimeStamp = graph.GetIntAttrDatE(edgeI, 'time')
			edgeTimeStamp = 1217567877
			currentEdgeSnapshotIndex = get_graph_index(edgeTimeStamp, currentEpochTime, timeQuantum)
			print currentEdgeSnapshotIndex
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

	"""
	Debug block

	g = graphSnapshots[3744]

	for e in g.Edges():
		print "e is " + str(e.GetSrcNId()) + " -> " + str(e.GetDstNId()) + ", time = " + str(g.GetIntAttrDatE(e.GetId(), 'time')) + ", edges = " + str(g.GetEdges())

	"""
	return graphSnapshots
	
def load_graph(version):
    # Loads the stack overflow full network without edge labels
    path = config.DATA_PATH[version]
    graph = snap.LoadEdgeList(snap.PNEANet, path, 0, 1, ' ')
    return graph

def print_graph_summary(graph):
    pass

def main():
    graph = load_graph('partial')
    print(graph.GetEdges())

    graph_snapshots = get_graph_snapshots(graph, 86400)

if __name__ == '__main__':
    main()
