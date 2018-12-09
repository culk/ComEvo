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
    # Load edge list and print sumary
    exp_name = '6months_12_leiden_exp'
    import_results = True
    graph = Graph('6months', degree=4)
    graph.print_summary()
    graph.set_num_time_slices(12)

    # Load or calculate communities and conductance
    if import_results:
        load_paths = {
                'communities': '../results/%s_communities.npy' % exp_name,
                'modularity': '../results/%s_modularity.npy' % exp_name,
                'sanitized_communities': '../results/%s_sanizited_communities.npy' % exp_name,
                'conductance': '../results/%s_conductance.npy' % exp_name,
            }
        graph.import_results(**load_paths)
    else:
        graph.calc_communities("leiden-algorithm", weight_fn=exp_fn, weighted=True)
        graph.sanitize_communities()
        graph.get_conductance(weight_fn=exp_fn)
        graph.export_results(exp_name)

    graph.algo_applied = exp_name

    # Plotting
    print('plotting modularity')
    graph.plot_modularity()
    print('getting best communities')
    best_comm_with_conductance = graph.select_best_communities(10)
    print('plotting conductance')
    graph.plot_conductance(best_comm_with_conductance)
    best_comm = map(lambda x: x[0], best_comm_with_conductance)
    comm_numnodes = graph.get_numnodes_from_comm_labels(best_comm)
    print('plotting size')
    graph.plot_numnodes(comm_numnodes)
    print('plotting clustering coefficient')
    graph.calc_clustering_coefficients()
    graph.plot_community_cc(best_comm)
    print('plotting egonet similarities')
    node_id = graph.select_best_egonet_node()
    graph.plot_egonet_community_similarity(node_id)
    print('plot egonets')
    graph.plot_egonet_community_similarity(node_id, distance=1)
    graph.plot_egonets(graph.algo_applied)

if __name__ == '__main__':
    main()
