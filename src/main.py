import pdb
import os
import snap
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import config

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

if __name__ == '__main__':
    main()
