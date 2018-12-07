import csv
import datetime
import time

import numpy as np
from scipy.cluster.vq import kmeans, whiten
import snap

import sys
import pdb

import config
import networkx as nx
from networkx import edge_betweenness_centrality
from networkx.algorithms import community

import igraph as ig
import leidenalg

class Graph():
    graph = None # snap.TNEANet containing all edges
    weights = None # List of weights for each time slice
    time_slices = [] # List of tuples containing [start, end) times
    sub_graphs = [] # List of snap.TUNGraph for each time slice
    communities = None # Matrix of shape T x N containing community labels --- For now considering as list of lists

    networkx_graph = None
    iGraph = None

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

        # Create a sorted edge list
        self.edge_list = []
        with open(path, 'r') as edge_list_file:
            for line in edge_list_file:
                src, dst, t = (int(i) for i in line.split(' '))
                # Ignore self edges
                if src == dst: continue
                self.edge_list.append((src, dst, t))
        self.edge_list = sorted(self.edge_list, key=lambda x: x[2])

        # Get start and end time stamps
        self._start_time = self.edge_list[0][2]
        self._end_time = self.edge_list[-1][2]

        # Initialize values for iterating over subgraphs
        self._cur_time = self._start_time
        self._cur_index = 0
        self.subgraph = snap.TNEANet.New()
        self.subgraph.AddFltAttrE('weight')

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

    def preprocess(self):
        snap.DelDegKNodes(self.graph, 1, 1)
        snap.DelDegKNodes(self.graph, 2, 2)
        snap.DelDegKNodes(self.graph, 3, 3)
        for edge in self.edge_list:
            # Parse the line
            src_id, dst_id, timestamp = edge

            # Add the nodes if not already present
            if not self.graph.IsNode(src_id): self.graph.AddNode(src_id)
            if not self.graph.IsNode(dst_id): self.graph.AddNode(dst_id)

            # Add the edge and assign the timestamp as an attribute
            edge_id = self.graph.AddEdge(src_id, dst_id)
            self.graph.AddIntAttrDatE(edge_id, timestamp, 'time')

        # Calculate self.node_to_index representing the whole self.graph
        self.node_to_index = dict()
        for i, n in enumerate(self.graph.Nodes()):
            self.node_to_index[n.GetId()] = i
    
    def calc_communities(self, method, weight_fn=None, weighted=False):
        """
        Update the internal weights, sub graphs, and communities by calling the
        function for community detection for the specified method.
        """
        # Calculate the bucket weights to use
        if weighted and weight_fn is not None:
            # TODO: should error here
            pass

        # Calculate community membership for each time slice
        if method == 'louvain':
            pass
        elif method == "lieden-algorithm":
            self.calc_communities_lieden_algorithm(weight_fn, weighted)
        elif method == 'girvan-newman':
            self.calc_communities_girvan_newman(weight_fn, weighted)
        elif method == 'spectral':
            k = 10 # TODO: replace with found best k
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

    def gen_next_subgraph(self, time_delta, weight_fn=None):
        start_time = self._cur_time
        end_time = start_time + time_delta
        i = self._cur_index
        while i < len(self.edge_list):
            src_id, dst_id, timestamp = self.edge_list[i]

            # Break if edge is outside of time slice
            if timestamp >= end_time:
                self._cur_index = i
                self._cur_time = end_time
                break

            # Add the nodes if not already present
            if not self.subgraph.IsNode(src_id): self.subgraph.AddNode(src_id)
            if not self.subgraph.IsNode(dst_id): self.subgraph.AddNode(dst_id)

            # Calculate the weight added by this edge
            if weight_fn:
                weight = weight_fn(self._start_time, self._end_time, timestamp)
            else:
                weight = 1

            # Add the edge weight to the subgraph
            if self.subgraph.IsEdge(src_id, dst_id):
                # Update the edge weight with the previous weight
                edge_id = self.subgraph.GetEI(src_id, dst_id).GetId()
                weight += self.subgraph.GetFltAttrDatE(edge_id, 'weight')
            else:
                # Add the edge and initialize with the weight
                edge_id = self.subgraph.AddEdge(src_id, dst_id)
            self.subgraph.AddFltAttrDatE(edge_id, weight, 'weight')
            i += 1

        # Calculate the time slice and return
        time_slice = (start_time, end_time)
        return self.subgraph, time_slice

    def update_subgraphs(self, time_delta):
        # TODO: debricate after switching to new subgraph function
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
            self.sub_graphs[i].AddIntAttrDatE(edge_id, timestamp, 'time')

    def save_subgraph_summaries(self, filename):
        # TODO: does not work with new subgraph function
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
        for i, sub_graph in enumerate(self.sub_graphs):

            #pdb.set_trace()

            networkxGraph = self.create_networkx_graph(sub_graph, weighted, weight_fn)

            self.networkx_graph = networkxGraph

            components = community.girvan_newman(networkxGraph, most_valuable_edge=self.most_central_edge)
            
            self.communities = tuple(sorted(component) for component in next(components))

            print str(self.communities)

            writeCommunityToFile(self.communities, i)

    def calc_communities_lieden_algorithm(self, weight_fn=None, weighted=False):
        """
        Create igraphs for each of the subgraphs so that the Lieden can work with it.
        Use Networkx Algorithm to Find Communities
        """
        #First do for the entire graph
        #pdb.set_trace()
        for graph_ind, sub_graph in enumerate(self.sub_graphs):

            #pdb.set_trace()

            iGraph = self.create_igraph(sub_graph, weighted, weight_fn)

            self.iGraph = iGraph

            #pdb.set_trace()

            to_delete_ids = [v.index for v in iGraph.vs if v.degree() == 0]
            iGraph.delete_vertices(to_delete_ids)

            #pdb.set_trace()

            partitions = leidenalg.find_partition(iGraph, leidenalg.ModularityVertexPartition, weights=iGraph.es["weight"]);
            
            #pdb.set_trace()

            print str(partitions)

            communityAssignment = {}
            for i in range(len(iGraph.vs)):
                if partitions.membership[i] in communityAssignment.keys():
                    communityAssignment[partitions.membership[i]].append(iGraph.vs[i]["name"])
                else:
                    communityAssignment[partitions.membership[i]] = [iGraph.vs[i]["name"]]

            self.writeParitionsToFile(communityAssignment, graph_ind)

    def calc_communities_spectral(self, k, weight_fn=None, weighted=False):
        for subgraph, time_slice in zip(self.sub_graphs, self.time_slices):
            print(time_slice)
            print(subgraph.GetNodes())
            print(subgraph.GetEdges())
            # TODO: implement weight function support
            # calculate normalized laplacian matrix and a node to index mapping
            D, A, node_to_index = self.get_matrices(subgraph, weight_fn)
            print(D)
            D_full = np.zeros_like(A)
            np.fill_diagonal(D_full, D)
            L = D_full - A
            D_star = np.zeros_like(A)
            np.fill_diagonal(D_star, D**(-.5))
            L_norm = np.dot(np.dot(D_star, L), D_star)

            # calculate top k eigenvectors
            values, vectors = np.linalg.eigh(L_norm)
            indices = np.argsort(values)[1:(k + 1)]
            features = vectors[:, indices]
            print(features)
            print(features.shape)

            # apply clustering algorithm
            features = whiten(features)
            centroids, distortion = kmeans(features, k)
            print(centroids)
            print(centroids.shape)
            print(distortion)

            # TODO: assign clusters based on closest centroid

    def get_matrices(self, graph, weight_fn=None):
        # TODO: update to use a weight function
        num_nodes = graph.GetNodes()
        D = np.zeros(num_nodes)
        A = np.zeros((num_nodes, num_nodes))

        node_to_index = dict()
        for i, n in enumerate(graph.Nodes()):
            node_to_index[n.GetId()] = i
            # degree is union of the set of out and in edges
            degree = len(set(x for x in n.GetOutEdges())
                         | set(y for y in n.GetInEdges()))
            D[i] = degree

        for e in graph.Edges():
            i = node_to_index[e.GetSrcNId()]
            j = node_to_index[e.GetDstNId()]
            A[i, j] = 1
            A[j, i] = 1

        return D, A, node_to_index

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

    def create_igraph(self, graph, weighted=False, weight_fn=None):
        """
        Creates a networkx graph for a given SNAP Graph
        """
        #pdb.set_trace()
        new_igraph = ig.Graph()
        edge_weights = []
        for edgeI in graph.Edges():
            edgeId = edgeI.GetId()
            srcId = edgeI.GetSrcNId()
            dstId = edgeI.GetDstNId()
            timestamp = graph.GetIntAttrDatE(edgeI, 'time')

            if srcId != dstId:
                #if not new_igraph.has_node(srcId):
                new_igraph.add_vertex(srcId)
                src_vertex_graph_index = new_igraph.vs.select(name=srcId)[0].index

                #if not networkxGraph.has_node(dstId):
                new_igraph.add_vertex(dstId)
                dst_vertex_graph_index = new_igraph.vs.select(name=dstId)[0].index

                #Add edges and weights
                #for now only keep 1 edge between nodes even if it has multi-edges
                if new_igraph.get_eid(src_vertex_graph_index, dst_vertex_graph_index, error=False) == -1:
                    if weighted:
                        edge_weight = weight_fn(self._start_time, self._end_time, timestamp)
                    else:
                        edge_weight = 1.0
                    new_igraph.add_edge(src_vertex_graph_index, dst_vertex_graph_index)
                    edge_weights.append(edge_weight)
        
        new_igraph.es['weight'] = edge_weights
        pdb.set_trace()
        return new_igraph

    def writeCommunityToFile(communities, index):
        communityAssignment = {}
        for i, community in enumerate(communities):
            communityAssignment[i] = sorted(community)

        with open("communityAssignment" + str(index) + ".txt", "w") as filename:
            filename.write(communityAssignment)
        filename.close()

    def writeParitionsToFile(self, communityAssignment, index):
        with open("communityAssignment" + str(index) + ".txt", "w") as filename:
            for key in communityAssignment.keys():
                filename.write(str(key) + ":" + str(communityAssignment[key]) + "\n\n")
        filename.close()
