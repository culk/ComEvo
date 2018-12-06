import csv
import datetime
import time

import numpy as np
import snap

import sys

import config
import networkx as nx
from networkx import edge_betweenness_centrality

class Graph():
    graph = None # snap.TNEANet containing all edges
    weights = None # List of weights for each time slice
    time_slices = [] # List of tuples containing [start, end) times
    sub_graphs = [] # List of snap.TUNGraph for each time slice
    communities = None # Matrix of shape T x N containing community labels --- For now considering as list of lists

    networkx_graph = None

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
    
    def calc_communities(self, method, weight_fn=None, weighted=False):
        """
        Update the internal weights, sub graphs, and communities by calling the
        function for community detection for the specified method.
        """
        # Calculate the bucket weights to use
        if weighted and weight_fn is not None:
            pass

        # Calculate community membership for each time slice
        if method == 'louvain':
            pass
        elif method == 'girvan-newman':
            self.calc_communities_girvan_newman(weight_fn, weighted)
        elif method == 'spectral':
            k = 100 # TODO: replace with found best k
            self.calc_communities_spectral(k, weight_fn, weighted)
        elif method == 'fastgreedy':
            self.calc_communities_fastgreedy()
        else: # etc.
            pass

    def get_conductance(self, weighted=False):
        """
        Calculates the average weighted or unweighted conductance for the graph
        given the already calculated sub graphs and communities.

        Returns a list of conductance values.
        """
        assert self.communities is not None

        minConductance = float(sys.maxint)
        graphNodes = list(self.networkx_graph.nodes)
        for community in self.communities:
            #conductance = nx.algorithms.cuts.conductance(self.networkx_graph, set(community), set(graphNodes).difference(set(community)), weight='weight')
            conductance = nx.algorithms.cuts.conductance(self.networkx_graph, set(community), set(graphNodes).difference(set(community)))
            print "conductance = " + str(conductance)
            if minConductance > conductance:
                minConductance = conductance
        return minConductance

    def print_summary(self):
        """
        Print summary statistics for the entire graph.
        """
        print('total nodes: %d' % self.graph.GetNodes())
        print('total edges: %d' % self.graph.GetEdges())
        print('first edge time: %d' % self._start_time)
        print('last edge time: %d' % self._end_time)
        print('time period: %d' % (self._end_time - self._start_time))

    def update_subgraphs(self, time_delta):
        """
        Update the self.sub_graphs to contain a list of graphs each containing
        only the edges present in that time slice.
        """
        # Initialize the sub graphs and time slices
        self.sub_graphs = []
        self.time_slices = []
        for t in range(self._start_time, self._end_time + 1, time_delta):
            self.sub_graphs.append(snap.TNEANet.New())
            self.time_slices.append((t, t + time_delta))

        # Add each edge to its respective sub graph
        for edge in self.graph.Edges():
            # Parse the edge
            src_id, dst_id = edge.GetSrcNId(), edge.GetDstNId()
            timestamp = self.graph.GetIntAttrDatE(edge, 'time')

            # Identify which subgraph index to add the edge to
            i = (timestamp - self._start_time) / time_delta

            # Add the nodes if not already present in the sub graph
            if not self.sub_graphs[i].IsNode(src_id):
                self.sub_graphs[i].AddNode(src_id)
            if not self.sub_graphs[i].IsNode(dst_id):
                self.sub_graphs[i].AddNode(dst_id)

            # Add the edge and assign the timestamp as an attribute, preserves
            # the edge id from the original graph.
            edge_id = self.sub_graphs[i].AddEdge(src_id, dst_id, edge.GetId())
            self.graph.AddIntAttrDatE(edge_id, timestamp, 'time')

    def save_subgraph_summaries(self, filename):
        """
        Write the summary statistics for each subgraph to a csv file.
        """
        field_names = [
                'slice_index',
                'start_time',
                'end_time',
                'num_nodes',
                'num_temporal_edges',
                'num_static_edges',
                ]

        # Open a csvfile and create a csv writer
        with open(filename, 'wb') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(field_names)

            # Calculate statistics for each sub graph
            for i, start_end in enumerate(self.time_slices):
                start_time, end_time = start_end
                field_values = [
                        i,
                        start_time,
                        end_time,
                        self.sub_graphs[i].GetNodes(),
                        self.sub_graphs[i].GetEdges(),
                        snap.CntUniqUndirEdges(self.sub_graphs[i]),
                        ]
                writer.writerow(field_values)

    #Private Mathods
    def calc_communities_girvan_newman(self, weight_fn=None, weighted=False):
        """
        Create Network graphs for each of the subgraphs
        Use Networkx Algorithm to Find Communities
        """
        #First do for the entire graph
        networkxGraph = self.create_networkx_graph(self.graph, weighted, weight_fn)

        self.networkx_graph = networkxGraph

        components = nx.algorithms.community.centrality.girvan_newman(networkxGraph, most_valuable_edge=self.most_central_edge)
        
        self.communities = tuple(sorted(component) for component in next(components))

        #print str(self.communities)

    def calc_communities_spectral(self, k, weight_fn=None, weighted=False):
        # for each subgraph (self.subgraph[i])
        # calculate laplacian matrix
        pass

    # TODO: Make sure community labels have consistent mapping across time slices
    def calc_communities_fastgreedy(self):
        communities = {}
        t_count = 0
        for subgraph in self.sub_graphs:
            subgraph_clean = snap.DelSelfEdges(subgraph)
            CmtyV = snap.TCnComV()
            modularity = snap.CommunityCNM(subgraph_clean, CmtyV)

            label_counter = 0
            for CnCom in CmtyV:
                for NI in CnCom:
                    nid = NI.GetId()
                    if nid not in communities:
                        communities[nid] = [0 for i in xrange(len(self.sub_graphs))]
                    communities[nid][t_count] = label_counter
                label_counter += 1
            t_count += 1
        self.communities = communities

    def most_central_edge(self, G):
        centrality = edge_betweenness_centrality(G, weight='weight')
        return max(centrality, key=centrality.get)

    def calculate_graph_modularity(self, graph, communities):
        assert communities is not None

        #Since communities is a list of lists, formulate snap vector and get the conductance
        communityModularities = []
        for community in communities:
            communityNodes = snap.TIntV()
            for nodeId in community:
                communityNodes.Add(nodeId)

            communityModularity = snap.GetModularity(graph, communityNodes)
            communityModularities.append(communityModularity)

        #Get the min of the communities - that is the modularity
        graphModularity = min(communityModularities)
        print str(graphModularity)
        return graphModularity

    def create_networkx_graph(self, graph, weighted=False, weight_fn=None):
        """
        Creates a networkx graph for a given SNAP Graph
        """
        networkxGraph = nx.Graph()
        for edgeI in graph.Edges():
            edgeId = edgeI.GetId()
            srcId = edgeI.GetSrcNId()
            dstId = edgeI.GetDstNId()
            timestamp = graph.GetIntAttrDatE(edgeI, 'time')

            if not networkxGraph.has_node(srcId):
                networkxGraph.add_node(srcId)

            if not networkxGraph.has_node(dstId):
                networkxGraph.add_node(dstId)

            #Add edges and weights
            #for now only keep 1 edge between nodes even if it has multi-edges
            if not networkxGraph.has_edge(srcId, dstId):
                if weighted:
                    edge_weight = weight_fn(self._start_time, self._end_time, timestamp)
                else:
                    edge_weight = 1.0
                networkxGraph.add_edge(srcId, dstId, weight=edge_weight)

        return networkxGraph
