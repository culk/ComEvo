import pdb
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
    # Select which version of the graph to load
    path = config.DATA_PATH[version]

    # Initialize a new graph and add an attribute for 'time'
    graph = snap.TNEANet.New()
    graph.AddIntAttrE('time')
    with open(path, 'r') as edge_list:
        for line in edge_list:
            # Parse the line
            src_id, dst_id, timestamp = (int(i) for i in line.split(' '))

            # Add the nodes if not already present
            if not graph.IsNode(src_id): graph.AddNode(src_id)
            if not graph.IsNode(dst_id): graph.AddNode(dst_id)

            # Add the edge and assign the timestamp as an attribute
            edge_id = graph.AddEdge(src_id, dst_id)
            graph.AddIntAttrDatE(edge_id, timestamp, 'time')

    return graph

def print_graph_summary(graph):
    print('total nodes: %d' % graph.GetNodes())
    print('total edges: %d' % graph.GetEdges())
    first_edge_time = graph.GetIntAttrDatE(0, 'time')
    last_edge_time = graph.GetIntAttrDatE(0, 'time')
    for edge in graph.Edges():
        edge_time = graph.GetIntAttrDatE(edge, 'time')
        if edge_time < first_edge_time:
            first_edge_time = edge_time
        if edge_time > last_edge_time:
            last_edge_time = edge_time
    print('first edge time: %d' % first_edge_time)
    print('last edge time: %d' % last_edge_time)
    print('time period: %d' % (last_edge_time - first_edge_time))

def main():
    graph = load_graph('partial')
    print_graph_summary(graph)
    pdb.set_trace()

    graph_snapshots = get_graph_snapshots(graph, 86400)

if __name__ == '__main__':
    main()
