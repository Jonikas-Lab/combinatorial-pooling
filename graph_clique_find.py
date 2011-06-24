#! /usr/bin/env python
"""
Clique-finding algorithm on graphs: given a graph, return a clique of size N if it exists.
Implements rudimentary Node and Graph classes (probably wouldn't be useful for a general purpose). 
 -- Weronika Patena, June 2011, Jonikas Lab, Carnegie Institution
"""

from itertools import combinations

class GraphNodeError(Exception):
    """ Exceptions in the clique_find module."""
    pass

class Node:
    """ A node has an arbitrary value, a set of neighbor nodes and a degree (equal to the number of neighbors)."""

    def __init__(self,value=None):
        """ Make a new node with the given value (default None), an empty neighbor set and a degree of 0."""
        self.value = value
        self.neighbors = set()
        self.degree = 0

    def add_neighbor(self,node):
        """ Add a neighbor node; raise GraphNodeError if it was already a neighbor."""
        if node in self.neighbors:  
            raise GraphNodeError("tried adding a neighbor that was already in this node's neighbor set!")
        self.degree += 1
        self.neighbors.add(node)

    def add_neighbor_nocheck(self,node):
        """ Add a neighbor node; do nothing if it was already a neighbor."""
        try:                add_neighbor(node)
        except GraphNodeError:  pass

    def remove_neighbor(self,node):
        """ Remove a neighbor node; raise GraphNodeError if it wasn't a neighbor."""
        if node not in self.neighbors:  
            raise GraphNodeError("tried removing a neighbor that wasn't in this node's neighbor set!")
        self.degree -= 1
        self.neighbors.remove(node)

    def remove_neighbor_nocheck(self,node):
        """ Remove a neighbor node; do nothing if it wasn't a neighbor."""
        try:                self.remove_neighbor(node)
        except GraphNodeError:  pass

    def clear_neighbors(self):
        """ Remove all neighbors of the node, and set its degree to 0."""
        self.neighbors.clear()
        self.degree = 0

    def check_degree(self):
        """ Raise GraphNodeError if self.degree isn't equal to the size of self.neighbors."""
        if not self.degree==len(self.neighbors):  raise GraphNodeError("node degree inconsistent!")

    def check_neighbors(self):
        """ Return True if the node is in the neighbor sets of all its neighbors, False otherwise.
        (Doesn't just raise an exception because sometimes this is a normal states for nodes removed from graphs)"""
        for node in self.neighbors:
            if self not in node.neighbors:  return False
        return True

    def __hash__(self):
        """ Use values for hashes - all values are assumed to be distinct, at least in one graph."""
        return hash(self.value)

    def __str__(self):
        """ A representation of a node is just the value."""
        return str(self.value)
    __repr__ = __str__ # TODO not sure if this is the best idea but works for now


class Graph:
    """ A graph is just a set of nodes and a size. All nodes are assumed to have distinct values."""

    def __init__(self,inputs=set()):
        """ Defines the nodes based on inputs (any sequence type): if an input isn't a node, it's used to make a node."""
        # TODO I know isinstance is bad but I'm not sure what else to sensibly do here
        self.nodes = set([x if isinstance(x,Node) else Node(x) for x in inputs])
        self.size = len(self.nodes)

    def build_edges_by_function(self,edge_function):
        """ A graph is just a set of nodes. Given a (node1,node2):if_edge function, fill in the neighbors for all nodes."""
        node_list = list(self.nodes)
        for (i,node1) in enumerate(node_list):
            for node2 in node_list[i+1:]:
                if edge_function(node1.value,node2.value):
                    node1.add_neighbor(node2)
                    node2.add_neighbor(node1)

    def add_edges_by_nodevalue(self,val_sets):
        """ For each set of values in val_sets, add an edge between the nodes with each pair of values in the set. 
        Probably extremely inefficient, for testing only. """
        val_to_node_dict = dict([(node.value,node) for node in self.nodes])
        for values in val_sets:
            for (val1,val2) in combinations(values,2):
                node1 = val_to_node_dict[val1]
                node2 = val_to_node_dict[val2]
                node1.add_neighbor(node2)
                node2.add_neighbor(node1)

    def clear(self):
        self.nodes = set()
        self.size = 0

    def remove_node(self,node1,keep_neighbor_record=False):
        """ Remove node1 from the graph node set and from the neighbor lists of all its neighbors; 
        remove all neighbors from its own neighbor list, unless keep_neighbor_record is set."""
        if node1 not in self.nodes: raise GraphNodeError("tried removing a node that wasn't in the graph!")
        for node2 in node1.neighbors:
            node2.remove_neighbor(node1)
        if not keep_neighbor_record:    
            node1.clear_neighbors()
        self.nodes.remove(node1)
        self.size -= 1

    def add_node_with_neighbors(self,node1):
        """ Add node1 to the graph node set and to the neighbor lists of all the nodes in its neighbor list."""
        if node1 in self.nodes: raise GraphNodeError("tried adding a node that was already in the graph!")
        for node2 in node1.neighbors:
            node2.add_neighbor(node1)
        self.nodes.add(node1)
        self.size += 1

    def remove_nodes_by_degree(self,min_degree,if_repeat=True):
        """ Remove nodes of degree below min_degree from graph. 
        Unless if_repeat is set to False, repeat until no such nodes are left."""
        N_removed = 0
        curr_nodes = self.nodes.copy()
        for node in curr_nodes:
            if node.degree<min_degree:   
                self.remove_node(node)
                N_removed += 1
        if if_repeat and N_removed:
            self.remove_nodes_by_degree(min_degree,if_repeat=True)

    def is_clique(self):
        """ Returns True if the graph is a clique, False otherwise."""
        max_degree = self.size-1
        for node in self.nodes:
            if not node.degree==max_degree:   return False
        return True

    def list_node_values(self):
        """ Return a set of the node values."""
        return set([n.value for n in self.nodes])

    def check_size(self):
        """ Raise GraphNodeError if self.size isn't equal to the size of self.nodes."""
        if not self.size==len(self.nodes):  raise GraphNodeError("graph size inconsistent!")

    def check_nodes(self):
        """ Raise GraphNodeError if any node has an incorrect degree or isn't in the neighbor sets of all its neighbors."""
        for node in self.nodes:     
            node.check_degree()
            if not node.check_neighbors():  raise GraphNodeError("graph node neighbor sets inconsistent!")

    def __str__(self):
        """ List of node values."""
        return "set(%s)"%', '.join([str(n) for n in self.nodes])
    __repr__ = __str__ # TODO not sure if this is the best idea but works for now

    def reduce_to_clique(self,min_size,required_nodes=None):
        """ Remove nodes from graph until it's a clique of size>=min_size; clear graph if there is no such clique.

        The final result is guaranteed to be either a clique of size at least min_size, or an empty graph; 
        it may be larger than min_size, but is not guaranteed to be the largest possible clique.
        The algorithm isn't brute-force, but it's unlikely to be optimal either.
        """
        if required_nodes is None: required_nodes = set()
        # remove all nodes with <N-1 edges (do this recursively until all nodes have at least N-1 edges)
        self.remove_nodes_by_degree(min_size-1,if_repeat=True)
        # if the graph is smaller than N, the task is impossible - clear the graph and fail
        if self.size<min_size:     
            self.clear()
            return
        # if the graph is a clique, we're done!
        if self.is_clique():   
            return
        # divide the nodes into ones that are definitely in the clique (max degree) and unknown.
        max_degree = self.size-1
        unknown_nodes = self.nodes - required_nodes
        required_nodes.update([node for node in unknown_nodes if node.degree==max_degree])
        unknown_nodes = self.nodes - required_nodes
        # If there's a clique, some nodes must not be in it. So try removing each node 
        #  and finding a clique in the sub-graph, recursively. 
        for node in unknown_nodes:
            # try removing one node, check if the smaller graph has an N-clique
            self.remove_node(node,keep_neighbor_record=True)
            self.reduce_to_clique(min_size,required_nodes)
            # if the result is a clique, we're done!
            if self.is_clique(): return
            # otherwise add the node back in and keep going (i.e. try removing a different node)
            self.add_node_with_neighbors(node)

# It would be nice to have a clique_find implementation that didn't modify the original graph...  NOT WORTH IT. 
#   I'd need to write custom deepcopy for my classes, and then copy the original but NOT at each recursion step... 
#   Anyway this is for one specific case of finding the clique, where I won't need the original back, so who cares.

# Thinking about implementation: if there is a clique, it must contain some nodes and not contain other nodes. so either start removing unknown_nodes from the graph and checking for cliqueness, or marking unknown_nodes as required_nodes (and then removing everything that's not its neighbor, and I suppose even removing it too and checking for a N-1 clique!) - try this for each node, and recurse down.  Which way would be better?  Depends on how the graph size compares to the clique size, really - if we have an N+1 graph and are looking for an N-clique, removing will be faster, but if we have a 4N graph, probably the other will be faster?  Anyway, removing seems cleaner, try that for now.

