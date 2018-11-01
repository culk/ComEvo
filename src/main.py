import os
import snap
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

DATA_DIR = '../data/'

def load_graph(version):
    # Loads the stack overflow full network without edge labels
    if version == 'full':
        path = os.path.join(DATA_DIR, 'stackoverflow_full.txt')
    elif version == 'partial':
        path = os.path.join(DATA_DIR, 'stackoverflow_partial.txt')
    graph = snap.LoadEdgeList(snap.PNEANet, path, 0, 1, ' ')
    return graph

def print_graph_summary(graph):
    pass

def main():
    graph = load_graph('partial')
    print(graph.GetEdges())

if __name__ == '__main__':
    main()
