from collections import Counter
import csv
import datetime
import pdb
import random
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
    modularity = None # Numpy array of shape T
    # Numpy matrix of shape C x T (where C is number of communities)
    conductance = None
    egonets = None
    egonet_node_id = None

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
        self.index_to_node = dict()
        for i, n in enumerate(self.graph.Nodes()):
            self.node_to_index[n.GetId()] = i
            self.index_to_node[i] = n.GetId()

    def sanitize_communities(self):
        """
        Works on the self.communities and sanities it to have the same community labels across timestamps
        """
        assert self.communities is not None

        sanitized = np.zeros_like(self.communities)
        sanitized.fill(-1)
        sanitized[:, 0] = self.communities[:, 0]
        last_label = 0

        for t in xrange(1, self.communities.shape[1]):
            prev_labels = set(sanitized[:, t-1])
            prev_labels.discard(-1)
            # Keep track of largest label assigned so far
            last_label = max(max(prev_labels), last_label)

            curr_labels = set(self.communities[:, t])
            curr_labels.discard(-1)

            # sort curr_labels by community size
            label_sizes = []
            for l in curr_labels:
                size = len(np.squeeze(np.where(self.communities[:, t] == l)))
                label_sizes.append((size, l))
            sorted_curr_labels = sorted(label_sizes, key=lambda x: -1 * x[0])
            
            for _, curr_l in sorted_curr_labels:
                #current timestep, get Jaccard wrt all the other communities
                curr_l_nodes = np.squeeze(np.where(self.communities[:, t] == curr_l))

                # Find label from previous time slice with greatest overlap
                max_score = 0.0
                new_l = -1
                for prev_l in prev_labels:
                    prev_l_nodes = np.squeeze(np.where(sanitized[:, t-1] == prev_l))

                    # Calculate jaccard similarity of these two communities
                    intersection = np.intersect1d(curr_l_nodes, prev_l_nodes)
                    union = np.union1d(curr_l_nodes, prev_l_nodes)
                    if len(union) != 0:
                        jaccard = len(intersection) / (1.0 * len(union))
                    else:
                        jaccard = 0.0
                    if jaccard > max_score:
                        max_score = jaccard
                        new_l = prev_l

                # If no overlapping community in previous time slice then assign
                # new labels.
                if new_l == -1:
                    last_label += 1
                    new_l = last_label

                # Update curr_l with new_l and remove new_l from prev_labels
                sanitized[curr_l_nodes, t] = new_l
                prev_labels.discard(new_l)

        self.sanitized_communities = sanitized
        self.communities = sanitized

        return sanitized

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
                    conductance[label, t] = sum_cut_edges[label] / float(denom)
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
            np.save('../results/%s_communities.npy' % exp_name, self.communities)
            print('Communities saved')
        if self.sanitized_communities is not None:
            print('Sanitized communities found, saving...')
            np.save('../results/%s_sanizited_communities.npy' % exp_name, self.sanitized_communities)
            print('Sanitized communities saved')
        if self.modularity is not None:
            print('Modularity found, saving...')
            np.save('../results/%s_modularity.npy' % exp_name, self.modularity)
            print('Modularity saved')
        if self.conductance is not None:
            print('Conductance found, saving...')
            np.save('../results/%s_conductance.npy' % exp_name, self.conductance)
            print('Conductance saved')

    def import_results(self, communities=None, sanitized_communities=None, modularity=None, conductance=None):
        if communities != None:
            self.communities = np.load(communities)
            print('Communities loaded')
        if sanitized_communities != None:
            self.sanitized_communities = np.load(sanitized_communities)
            print('Sanitized communities loaded')
        if modularity != None:
            self.modularity = np.load(modularity)
            print('Modularity loaded')
        if conductance != None:
            self.conductance = np.load(conductance)
            print('Conductance loaded')

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

        np.save('fastgreedy-assignments.npy', communities)
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

    def select_best_communities(self, num_best=10):
        conductance = zip(range(len(self.conductance)), self.conductance)
        # Filter out all-zero rows
        conductance = filter(lambda x: len(np.nonzero(x[1])[0]) > 0, conductance)
        # Filter out rows with last timeslice conductance == 0
        final_nonzero_conductance = filter(lambda x: x[1][-1] > 0, conductance)
        best_conductance = sorted(final_nonzero_conductance, key=lambda x: x[1][-1], reverse=False)[:num_best]
        # Sort best_conductance again based on community label ID
        best_conductance = sorted(best_conductance, key=lambda x: x[0])
        return best_conductance

    def get_numnodes_from_comm_labels(self, comm_labels):
        counts_at_time = []
        for i in xrange(self.num_time_slices):
            unique, counts = np.unique(self.communities[:, i], return_counts=True)
            comm_counts = dict(zip(unique, counts))
            counts_at_time.append(comm_counts)
        comm_numnodes = []
        for comm_label in comm_labels:
            comm_label_size = np.zeros((self.num_time_slices))
            for i in xrange(self.num_time_slices):
                if comm_label in counts_at_time[i]:
                    comm_label_size[i] = counts_at_time[i][comm_label]
            comm_numnodes.append([comm_label, comm_label_size])
        return comm_numnodes

    def plot_modularity(self):
        plt.figure()
        plt.plot(range(len(self.modularity)), self.modularity, 'o-')
        plt.xlabel('Cumulative Time Slice #')
        plt.ylabel('Graph Modularity')
        plt.title('Temporal Community Evolution - Graph Modularity (%s)' % self.algo_applied)
        plt.savefig('modularity_%s.png' % self.algo_applied)

    def plot_conductance(self, y_data=None, max_communities=10):
        """
        y_data: list[] of 2-tuples
            (community_label, np.array(conductance_values_at_timeslices))
        """
        plot_colors = [
            'blue', 'green', 'red', 'cyan', 'magenta',
            'grey', 'black', 'pink', 'purple', 'orange'
        ]
        if y_data == None:
            y_data = zip(range(len(self.conductance)), self.conductance)
        plt.figure()
        cond_plot = plt.subplot(111)
        fontP = FontProperties()
        fontP.set_size('small')
        for i in xrange(min(max_communities, len(y_data))):
            cond_plot.plot(range(len(y_data[i][1])), y_data[i][1], 'o-', label=y_data[i][0], color=plot_colors[i % len(plot_colors)])
        cond_plot.set_xlabel('Cumulative Time Slice #')
        cond_plot.set_ylabel('Community Conductance')
        cond_plot.set_title('Temporal Community Evolution - Conductance (%s)' % self.algo_applied)
        cond_plot.legend(loc="upper left", bbox_to_anchor=(1,1), prop=fontP)
        plt.savefig('conductance_%s.png' % self.algo_applied)

    def select_best_egonet_node(self):
        assert self.conductance is not None

        # Find a label for a community with the best final conductance that
        # existed in the largest number of time slices.
        timeslice = 0
        label = -1
        best_conductance = 1
        for i in xrange(self.communities.shape[1] - 1):
            communities = set(self.communities[:, i]) | set(self.communities[:, -1])
            if len(communities) > 0:
                timeslice = i
                break
        for l in communities:
            curr_conductance = self.conductance[l, -1]
            if curr_conductance != 0 and curr_conductance < best_conductance:
                best_conductance = curr_conductance
                label = l

        # Count node membership across all communities
        node_membership_count = Counter()
        for t in xrange(timeslice, self.communities.shape[1]):
            for n in np.squeeze(np.where(self.communities[:, t] == label)):
                node_id = self.index_to_node[n]
                node_membership_count[node_id] += 1

        # Return the node with the best membership
        return node_membership_count.most_common(1)[0][0]

    def plot_egonet_community_similarity(self, node_id, distance=2):
        assert self.communities is not None
        assert node_id != -1

        community_similarity = []
        self.egonets = []
        self.egonet_edge_lists = []
        self.egonet_node_id = node_id

        # Calculate the similarity scores for each subgraph
        t = 0
        for subgraph, _ in self.gen_next_subgraph():
            # Calculate the egonet set
            egonet = set([node_id]) # visited
            node = subgraph.GetNI(node_id)
            curr_set = set([node]) # to visit
            for _ in xrange(distance):
                next_set = set()
                for n in curr_set:
                    # Add all 'in' neighbors to the egonet
                    for i in xrange(n.GetInDeg()):
                        neighbor_id = n.GetInNId(i)
                        if neighbor_id not in egonet:
                            neighbor = subgraph.GetNI(neighbor_id)
                            next_set.add(neighbor)
                            egonet.add(neighbor_id)
                    # Add all 'out' neighbors to the egonet
                    for i in xrange(n.GetOutDeg()):
                        neighbor_id = n.GetOutNId(i)
                        if neighbor_id not in egonet:
                            neighbor = subgraph.GetNI(neighbor_id)
                            next_set.add(neighbor)
                            egonet.add(neighbor_id)
                curr_set = next_set

            # Append the percent similar for current time slice to scores
            num_similar = 0
            for n in egonet:
                a = self.node_to_index[n]
                b = self.node_to_index[node_id]
                if self.communities[a, t] == self.communities[b, t]:
                    num_similar += 1
            similarity = float(num_similar) / float(len(egonet))
            community_similarity.append(similarity)

            # Calculate the list of egonet edges
            egonet_edge_list = []
            for src in egonet:
                for dst in egonet:
                    if subgraph.IsEdge(src, dst):
                        egonet_edge_list.append((src, dst))

            self.egonets.append(egonet)
            self.egonet_edge_lists.append(egonet_edge_list)
            t += 1

        # Plot the values
        plt.figure()
        x = range(len(community_similarity))
        y = community_similarity
        plt.plot(x, y, 'o-')
        for t, xy in enumerate(zip(x, y)):
            l = self.communities[self.node_to_index[node_id], t]
            plt.annotate(l, xy=xy)
        plt.xlabel('Cumulative Time Slice #')
        plt.ylabel('Same community (%)')
        plt.title("Percent nodes in Node %d's Egonet with Same Community (distance = %d)" % (node_id, distance))
        plt.savefig('community_similarity_node%d_%s.png' % (node_id, self.algo_applied))

    def calc_clustering_coefficients(self):
        assert self.communities is not None

        # Initialize clustering coefficien values.
        num_communities = np.max(self.communities) + 1
        clustering_coefficients = np.zeros((num_communities, self.num_time_slices))

        # Calculate the clustering coefficiens for each time slice.
        t = 0
        for subgraph, time_slice in self.gen_next_subgraph():
            for l in xrange(num_communities):
                node_ccs = []
                nodes = np.atleast_1d(np.squeeze(np.where(self.communities[:, t] == l)))
                for node_index in nodes:
                    node_id = self.index_to_node[node_index]
                    node_ccs.append(snap.GetNodeClustCf(subgraph, node_id))

                if len(node_ccs) == 0:
                    community_cc = 0.0
                else:
                    community_cc = sum(node_ccs) / float(len(node_ccs))
                clustering_coefficients[l, t] = community_cc
            t += 1

        self.clustering_coefficients = clustering_coefficients
        return clustering_coefficients

    def plot_community_cc(self, communities=None, max_communities=10):
        """
        y_data: list[] of 2-tuples
            (community_label, np.array(clustering_coefficient_at_timeslices))
        """
        plot_colors = [
            'blue', 'green', 'red', 'cyan', 'magenta',
            'grey', 'black', 'pink', 'purple', 'orange'
        ]

        if communities == None:
            y_data = zip(range(len(self.clustering_coefficients)), self.clustering_coefficients)
        else:
            y_data = zip(communities, self.clustering_coefficients[communities])
        plt.figure()
        cond_plot = plt.subplot(111)
        fontP = FontProperties()
        fontP.set_size('small')
        for i in xrange(min(max_communities, len(y_data))):
            cond_plot.plot(range(len(y_data[i][1])), y_data[i][1], 'o-', label=y_data[i][0], color=plot_colors[i % len(plot_colors)])
        cond_plot.set_xlabel('Cumulative Time Slice #')
        cond_plot.set_ylabel('Clustering Coefficient')
        cond_plot.set_title('Temporal Community Evolution - Clustering Coefficient (%s)' % self.algo_applied)
        cond_plot.legend(loc="upper left", bbox_to_anchor=(1,1), prop=fontP)
        plt.savefig('clustering_coefficient_%s.png' % self.algo_applied)

    def plot_numnodes(self, comm_numnodes):
        """
        y_data: list[] of 2-tuples
            (community_label, np.array(conductance_values_at_timeslices))
        """
        plot_colors = [
            'blue', 'green', 'red', 'cyan', 'magenta',
            'grey', 'black', 'pink', 'purple', 'orange'
        ]
        plt.figure()
        numnodes_plot = plt.subplot(111)
        fontP = FontProperties()
        fontP.set_size('small')
        for i in xrange(len(comm_numnodes)):
            numnodes_plot.semilogy(range(len(comm_numnodes[i][1])), comm_numnodes[i][1], 'o-', label=comm_numnodes[i][0], color=plot_colors[i % len(plot_colors)])
        numnodes_plot.set_xlabel('Cumulative Time Slice #')
        numnodes_plot.set_ylabel('Community Size')
        numnodes_plot.set_title('Temporal Community Evolution - Community Size (%s)' % self.algo_applied)
        numnodes_plot.legend(loc="upper left", bbox_to_anchor=(1,1), prop=fontP)
        plt.savefig('numnodes_%s.png' % self.algo_applied)
