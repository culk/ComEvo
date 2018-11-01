import datetime
import time

import numpy as np
import snap

import config


class Graph():
    graph = None # snap.TNEANet containing all edges
    weights = None # List of weights for each time slice
    sub_graphs = [] # List of snap.TUNGraph for each time slice
    communities = None # Matrix of shape T x N containing community labels

    # Unix timestamp of first and last edge
    _start_time = None
    _end_time = None

    def __init__(self, version):
        """
        Takes in a string representing which graph version to use and loads
        the entire edge list into a snap.TNEANet.
        """
        # Select which version of the graph to load
        path = config.DATA_PATH[version]

        # Initialize a new graph and add an attribute for 'time'
        self.graph = snap.TNEANet.New()
        self.graph.AddIntAttrE('time')
        with open(path, 'r') as edge_list:
            for line in edge_list:
                # Parse the line
                src_id, dst_id, timestamp = (int(i) for i in line.split(' '))

                # Update start and end time
                if not self._start_time or self._start_time > timestamp:
                    self._start_time = timestamp
                if not self._end_time or self._end_time < timestamp:
                    self._end_time = timestamp

                # Add the nodes if not already present
                if not self.graph.IsNode(src_id): self.graph.AddNode(src_id)
                if not self.graph.IsNode(dst_id): self.graph.AddNode(dst_id)

                # Add the edge and assign the timestamp as an attribute
                edge_id = self.graph.AddEdge(src_id, dst_id)
                self.graph.AddIntAttrDatE(edge_id, timestamp, 'time')
    
    def calc_communities(self, method, time_delta, weight_fn=None, weighted=False):
        """
        Update the internal weights, subgraphs, and communities by calling the
        function for community detection for the specified method.
        """
        # Calculate the bucket weights to use
        if weighted and weight_fn is not None:
            pass

        # Organize the subgraphs by the time_delta

        # Calculate community membership for each time slice
        if method == 'louvain':
            pass
        else: # etc.
            pass

    def get_conductance(self, weighted=False):
        """
        Calculates the average weighted or unweighted conductance for the graph
        given the already calculated subgraphs and communities.

        Returns a list of conductance values.
        """
        assert communities is not None

        pass

    def print_summary(self):
        """
        Print summary statistics for the entire graph.
        """
        print('total nodes: %d' % self.graph.GetNodes())
        print('total edges: %d' % self.graph.GetEdges())
        print('first edge time: %d' % self._start_time)
        print('last edge time: %d' % self._end_time)
        print('time period: %d' % (self._end_time - self._start_time))

