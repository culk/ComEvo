from collections import Counter
import csv
import datetime
import pdb
import sys
import time

import igraph as ig
import leidenalg
import networkx as nx # TODO: remove when conductance updated
import numpy as np
import snap

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties

import config


class Graph():
    graph = None # snap.TNEANet containing all edges
    sanitized_communities = None # Matrix of shape N * T containing community labels
    communities = None # Numpy matrix of shape N x T containing community labels
    sanitized_communities = None # Numpy matrix of shape N x T containing community labels
    modularity = None # Numpy array of shape T
    # Numpy matrix of shape C x T (where C is number of communities)
    conductance = None

    edge_list = []

    networkx_graph = None
    iGraph = None

    # Unix timestamp of first and last edge
    _start_time = None
    _end_time = None

    # Attributes for subgraphs
    time_delta = None
    subgraph = None # snap.TNEANet containing edges up to current time slice
    _cur_time = None
    _cur_index = 0

    def __init__(self, version, degree=0):
        """
        Takes in a string representing which graph version to use and loads
        the entire edge list into a snap.TNEANet.

            version: version of graph edge list to use (see config).
            degree: remove nodes from edge list with total degree less than this.
        """
        # Select which version of the graph to load
        path = config.DATA_PATH[version]

        # Create a sorted edge list, removing low degree nodes
        self.edge_list = []
        degree_counts = Counter()
        with open(path, 'r') as edge_list_file:
            for line in edge_list_file:
                src, dst, t = (int(i) for i in line.split(' '))
                # Ignore self edges
                if src == dst: continue
                degree_counts[src] += 1
                degree_counts[dst] += 1
                self.edge_list.append((src, dst, t))
        if degree > 0:
            to_remove = [node_id for node_id, d in degree_counts.iteritems() if d <= degree]
            new_edge_list = []
            for edge in self.edge_list:
                src, dst, _ = edge
                if src in to_remove or dst in to_remove:
                    continue
                else:
                    new_edge_list.append(edge)
            self.edge_list = new_edge_list
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

        # Indicator of which community detection algorithm has been applied
        self.algo_applied = 'None'

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

    def sanitize_communities(self):
        """
        Works on the self.communities and sanities it to have the same community labels across timestamps
        """
        assert self.communities is not None

        sanitized_communities = np.zeros_like(self.communities)
        sanitized_communities.fill(-2)

        num_nodes = np.shape(self.communities)[0]
        num_timesteps = np.shape(self.communities)[1]

        #set the first one as is

        global_counter = 0

        sanitized_communities[:, 0] = self.communities[:, 0]

        for timestep in range(0, num_timesteps-1):
            #go through each community at the current timestep and try to sanitize the next timestep based on the Jaccard Index match
            current_distinct_communities = np.unique(sanitized_communities[:, timestep])
            current_distinct_communities = np.delete(current_distinct_communities, np.where(current_distinct_communities == -2))
            current_distinct_communities = np.delete(current_distinct_communities, np.where(current_distinct_communities == -1))

            next_distinct_communities = np.unique(self.communities[:, timestep+1])
            next_distinct_communities = np.delete(next_distinct_communities, np.where(next_distinct_communities == -1))
            
            for i in range(len(current_distinct_communities)):
                #current timestep, get Jaccard wrt all the other communities
                nodes_i = np.where(sanitized_communities[:, timestep] == current_distinct_communities[i])[0]

                scores = []
                for j in range(len(next_distinct_communities)):
                    #get elements belonging to communities j (in next timestep)
                    nodes_j = np.where(self.communities[:, timestep+1] == next_distinct_communities[j])[0]

                    intersection = np.intersect1d(nodes_i, nodes_j)
                    union = np.union1d(nodes_i, nodes_j)

                    if len(union) != 0:
                        jaccard = len(intersection) / (1.0 * len(union))
                    else:
                        jaccard = 0.0

                    scores.append(jaccard)

                #if the elements in next_distinct_communities reached 0, we are already done
                if len(next_distinct_communities) != 0:
                    next_label_i = next_distinct_communities[np.argmax(scores)]

                    #assign elements to max jaccard
                    indices_to_update = np.argwhere(self.communities[:, timestep + 1] == next_label_i)

                    sanitized_communities[np.squeeze(indices_to_update), timestep + 1] = current_distinct_communities[i]

                    #delete the community label to not process any further
                    next_distinct_communities = np.delete(next_distinct_communities, np.argwhere(next_distinct_communities == next_label_i))

            global_counter = max(i, global_counter + 1)

            #if there are more elements in j, then update them
            for j in range(len(next_distinct_communities)):
                indices_to_update = np.argwhere(self.communities[:, timestep + 1] == next_distinct_communities[j])
                sanitized_communities[np.squeeze(indices_to_update), timestep + 1] = global_counter
                global_counter += 1

        sanitized_communities[sanitized_communities == -2] = -1
        self.sanitized_communities = sanitized_communities

        return sanitized_communities

    def calc_communities(self, method, weight_fn=None, weighted=False):
        """
        Update the internal weights, sub graphs, and communities by calling the
        function for community detection for the specified method.
        """
        # Calculate community membership for each time slice
        if method == "leiden-algorithm":
            self.calc_communities_leiden_algorithm(weight_fn)
        elif method == 'fastgreedy':
            self.calc_communities_fastgreedy()
        else: # etc.
            pass
        self.algo_applied = method

    def get_conductance(self, weight_fn=None):
        """
        Calculates the weighted or unweighted conductance for the graph for each
        community at each time slice given the already calculated communities
        for each subgraph.

        Returns a matrix of conductance values.
        """
        assert self.communities is not None

        # Initialize conductance values.
        num_communities = np.max(self.communities) + 1
        conductance = np.zeros((num_communities, self.num_time_slices))

        # Calculate the conductance for each time slice.
        t = 0
        for subgraph, time_slice in self.gen_next_subgraph(weight_fn):
            # Track the sum of the edge weights for calculating conductance for
            # each community. Communities that do not exist in the time slice
            # will have conductance of 0.
            volume_S_edges = np.zeros(num_communities)
            sum_cut_edges = np.zeros(num_communities)
            for edge in subgraph.Edges():
                # Calculate the community labels for each node.
                src_id, dst_id = edge.GetSrcNId(), edge.GetDstNId()
                weight = subgraph.GetFltAttrDatE(edge, 'weight')
                src_community = self.communities[self.node_to_index[src_id], t]
                dst_community = self.communities[self.node_to_index[dst_id], t]

                # All nodes in the current subgraph should have a community.
                assert src_community != -1 and dst_community != -1

                # For calculating the numerator (the edges originating from the
                # community).
                if src_community != dst_community:
                    sum_cut_edges[src_community] += weight
                # For calculating the denominator (min of the sum of edges
                # originating from the community or from the rest of the graph).
                volume_S_edges[src_community] += weight

            # Use above values for calculating conductance for each community.
            for label in xrange(num_communities):
                denom = min(volume_S_edges[label], np.sum(volume_S_edges) - volume_S_edges[label])
                if sum_cut_edges[label] == 0:
                    conductance[label, t] = 0
                else:
                    conductance[label, t] = sum_cut_edges[label] / denom
            t += 1

        self.conductance = conductance

        return conductance

    def print_summary(self):
        """
        Print summary statistics for the entire graph.
        """
        print('total nodes: %d' % self.graph.GetNodes())
        print('total edges: %d' % self.graph.GetEdges())
        print('first edge time: %d' % self._start_time)
        print('last edge time: %d' % self._end_time)
        print('time period: %d' % (self._end_time - self._start_time))

    def export_results(self, exp_name='test'):
        """
        Save a copy of the values calculated for the current experiment.

        exp_name should be a unique string used to identify the files saved.
        """
        if self.communities is not None:
            print('Communities found, saving...')
            np.save('../results/%s_communities.npy', self.communities)
            print('Communities saved')
        if self.sanitized_communities is not None:
            print('Sanitized communities found, saving...')
            np.save('../results/%s_sanizited_communities.npy', self.sanitized_communities)
            print('Sanitized communities saved')
        if self.modularity is not None:
            print('Modularity found, saving...')
            np.save('../results/%s_modularity.npy', self.modularity)
            print('Modularity saved')
        if self.conductance is not None:
            print('Conductance found, saving...')
            np.save('../results/%s_conductance.npy', self.conductance)
            print('Conductance saved')

    def set_time_delta(self, time_delta):
        """
        Update the time delta and time slices using the given time_delta.
        """
        self.time_delta = time_delta
        self.num_time_slices = (self._end_time - self._start_time + 1) / time_delta

    def set_num_time_slices(self, num_time_slices):
        """
        Update the time delta and time slices using the given num_time_slices.
        """
        self.num_time_slices = num_time_slices
        self.time_delta = (self._end_time - self._start_time + 1) / num_time_slices

    def gen_next_subgraph(self, weight_fn=None):
        """
        Generate pairs of (subgraph, time slice) which contains snapshots of
        the graph based on the edges known at the end of the current time slice.
        """
        # Reset subgraph state
        self._cur_index = 0
        self._cur_time = self._start_time
        self.subgraph.Clr()
        del self.subgraph
        self.subgraph = snap.TNEANet.New()
        self.subgraph.AddFltAttrE('weight')

        # Initialize time values
        start_time = self._cur_time
        end_time = start_time + self.time_delta
        i = self._cur_index
        while i < len(self.edge_list):
            # Parse the edge
            src_id, dst_id, timestamp = self.edge_list[i]

            # Yield if edge is outside of time slice
            if timestamp >= end_time:
                time_slice = (start_time, end_time)
                self._cur_index = i
                self._cur_time = end_time
                start_time = self._cur_time
                end_time = start_time + self.time_delta
                yield self.subgraph, time_slice

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
                if weight_fn is not None:
                    weight += self.subgraph.GetFltAttrDatE(edge_id, 'weight')
            else:
                # Add the edge and initialize with the weight
                edge_id = self.subgraph.AddEdge(src_id, dst_id)
            self.subgraph.AddFltAttrDatE(edge_id, weight, 'weight')
            i += 1

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
            i = 0
            for subgraph, time_slice in self.gen_next_subgraph():
                start_time, end_time = time_slice
                field_values = [
                        i,
                        start_time,
                        end_time,
                        subgraph.GetNodes(),
                        subgraph.GetEdges(),
                        ]
                writer.writerow(field_values)
                i += 1

    #Private Mathods
    def calc_communities_leiden_algorithm(self, weight_fn=None):
        """
        Create igraphs for each of the subgraphs so that the Lieden can work with it.
        Use Networkx Algorithm to Find Communities
        """
        # Initialize communities
        communities = -1 * np.ones((len(self.node_to_index), self.num_time_slices), dtype=int)
        modularity = np.zeros((self.num_time_slices), dtype=float)

        t = 0
        for subgraph, time_slice in self.gen_next_subgraph(weight_fn):
            # Create an igraph representation of the subgraph
            self.iGraph = self.create_igraph(subgraph)

            # Delete degree 0 nodes
            to_delete_ids = [v.index for v in self.iGraph.vs if v.degree() == 0]
            self.iGraph.delete_vertices(to_delete_ids)

            # Calculate communities
            partitions = leidenalg.find_partition(self.iGraph,leidenalg.ModularityVertexPartition, weights=self.iGraph.es['weight'])
            modularity[t] = partitions.quality()
            for i in xrange(len(self.iGraph.vs)):
                node_id = self.iGraph.vs[i]['name']
                communities[self.node_to_index[node_id], t] = partitions.membership[i]

            print('Time slice: %s, Modularity = %f' % (time_slice, modularity[t]))
            t += 1

        np.save('leiden-assignments.npy', communities)
        self.communities = communities
        self.modularity = modularity

    def calc_communities_fastgreedy(self):
        """
        Use snap's implementation of fastgreedy algorithmn to calculate
        community labels for each node in every time slice subgraph.
        """
        communities = -1 * np.ones((len(self.node_to_index), self.num_time_slices), dtype=int)
        modularity = np.zeros((self.num_time_slices), dtype=float)

        t = 0
        for subgraph, time_slice in self.gen_next_subgraph():
            subgraph_clean = snap.ConvertGraph(snap.PUNGraph, subgraph)
            snap.DelSelfEdges(subgraph_clean)
            CmtyV = snap.TCnComV()
            modularity[t] = snap.CommunityCNM(subgraph_clean, CmtyV)

            for label, CnCom in enumerate(CmtyV):
                for node_id in CnCom:
                    communities[self.node_to_index[node_id], t] = label

            print('Time slice: %s, Modularity = %f' % (time_slice, modularity[t]))
            t += 1

        #for i in xrange(len(communities)):
            #if communities[i, :3] != [-1, -1, -1]:
                #print communities[i][:3]

        self.communities = communities
        self.modularity = modularity

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

    def create_igraph(self, graph):
        """
        Creates a networkx graph for a given SNAP Graph
        """
        new_igraph = ig.Graph()
        edge_weights = []
        for edgeI in graph.Edges():
            # Parse the edge information 
            edgeId = edgeI.GetId()
            srcId = edgeI.GetSrcNId()
            dstId = edgeI.GetDstNId()
            weight = graph.GetFltAttrDatE(edgeI, 'weight')

            # Add the nodes if not already present
            if len(new_igraph.vs) == 0 or srcId not in new_igraph.vs['name']:
                new_igraph.add_vertex(srcId)
            if len(new_igraph.vs) == 0 or dstId not in new_igraph.vs['name']:
                new_igraph.add_vertex(dstId)

            # Add edge and weight
            src_vertex_graph_index = new_igraph.vs.select(name=srcId)[0].index
            dst_vertex_graph_index = new_igraph.vs.select(name=dstId)[0].index
            new_igraph.add_edge(src_vertex_graph_index, dst_vertex_graph_index)
            edge_weights.append(weight)
        
        new_igraph.es['weight'] = edge_weights
        return new_igraph

    def writeCommunityToFile(self, communities, index):
        communityAssignment = {}
        for i, community in enumerate(communities):
            communityAssignment[i] = sorted(community)

        with open("communityAssignment" + str(index) + ".txt", "w") as filename:
            filename.write(communityAssignment)
        filename.close()

    def writeParitionsToFile(self, communityAssignment, index):
        with open("communityAssignment" + str(index) + ".txt", "w") as filename:
            for key in communityAssignment.keys():
                filename.write(str(key) + ":" + str(communityAssignment[key]) + "\n")
        filename.close()

    def plot_modularity(self):
        plt.plot(range(1, len(self.modularity)+1), self.modularity, 'o-')
        plt.xlabel('Cumulative Time Slice #')
        plt.ylabel('Grpah Modularity')
        plt.title('Temporal Community Evolution - Graph Modularity (%s)' % self.algo_applied)
        plt.savefig('modularity_%s.png' % self.algo_applied)

    def plot_conductance(self, y_data=self.conductance, max_communities=10):
        cond_plot = plt.subplot(111)
        fontP = FontProperties()
        fontP.set_size('small')
        for i in xrange(min(max_communities, len(y_data))):
            cond_plot.plot(range(1, len(y_data[i])+1), y_data[i], 'o-', label=str(i))
        cond_plot.set_xlabel('Cumulative Time Slice #')
        cond_plot.set_ylabel('Community Conductance')
        cond_plot.set_title('Temporal Community Evolution - Conductance (%s)' % self.algo_applied)
        cond_plot.legend(loc="upper left", bbox_to_anchor=(1,1), prop=fontP)
        plt.savefig('conductance_%s.png' % self.algo_applied)
