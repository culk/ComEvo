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
    graph = Graph('1year')
    #graph.print_summary()
    #graph.preprocess() # only needed to delete small nodes
    graph.set_num_time_slices(12)
    graph.calc_communities("leiden-algorithm", weight_fn=exp_fn, weighted=True)
    #graph.calc_communities('fastgreedy', weight_fn=None, weighted=False)
    graph.export_results('leiden_base')
    graph.sanitize_communities()
    graph.get_conductance(weight_fn=exp_fn)
    graph.export_results('leiden_sanitized')
    #graph.plot_modularity()

if __name__ == '__main__':
    main()
