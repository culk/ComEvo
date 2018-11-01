import os
import snap
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

DATA_DIR = 'data/'

def load_graph():
    # Loads the stack overflow full network without edge labels
    path = os.path.join(DATA_DIR, 'stackoverflow_full.txt')
    graph = snap.LoadEdgeList(snap.PNEANet, path, 0, 1, ' ')

def print_graph_summary(graph):
    pass

def main():
    graph = load_graph()


if __name__ == '__main__':
    main()
