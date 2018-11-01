import os
import snap
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import config

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

if __name__ == '__main__':
    main()
