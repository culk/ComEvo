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
    graph = Graph('30days')
    graph.print_summary()
    graph.preprocess()
    graph.set_num_time_slices(10)
    graph.calc_communities("leiden-algorithm", weight_fn=linear_fn, weighted=True)
    #graph.calc_communities("fastgreedy", weight_fn=None, weighted=False)
    #graph.get_conductance()

if __name__ == '__main__':
    main()
