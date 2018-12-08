import datetime
import os
import pdb
import sys
import time

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import snap

import config
from graph import Graph
from utilities import *


def main():
    graph = Graph('6months', degree=4)
    graph.print_summary()
    graph.set_num_time_slices(12)
    graph.calc_communities("leiden-algorithm", weight_fn=exp_fn, weighted=True)

    '''
    load_paths = [
        '../results/6months_leiden_communities.npy',
        '../results/6months_leiden_sanizited_communities.npy',
        '../results/6months_leiden_modularity.npy',
        '../results/6months_leiden_conductance.npy'
    ]
    graph.import_results(*load_paths)
    '''

    #graph.calc_communities('fastgreedy', weight_fn=None, weighted=False)
    graph.export_results('6months_leiden')
    graph.sanitize_communities()
    graph.export_results('6months_leiden')
    graph.get_conductance(weight_fn=exp_fn)
    graph.export_results('6months_leiden')
    #graph.plot_modularity()
    #best_comm_with_conductance = graph.select_best_communities(10)
    #graph.plot_conductance(best_comm_with_conductance)
    #best_comm = map(lambda x: x[0], best_comm_with_conductance)
    #comm_numnodes = graph.get_numnodes_from_comm_labels(best_comm)
    #graph.plot_numnodes(comm_numnodes)

if __name__ == '__main__':
    main()
