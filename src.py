#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 12:41:14 2023

@author: tezcan.b
"""

import timeit 
import numpy as np
from itertools import product, repeat
import networkx as nx 
import gurobipy as gp
from gurobipy import GRB
import pandas as pd 

#Class/Object of our problem instance 

class Network:
      #reads capacities, interdiction probabilities and cost for the specified network
    #note that network id and data location should match 
    def __init__(self, n_id : 'str',
                 folder = '/Data'):
        
        import networkx as nx 
        import pandas as pd
        import json 
        
        self.n_id = n_id
        self.folder = folder
        self.cap = pd.read_fwf(folder+'/'+n_id+'/'+'cap.txt', header = None)
        self.prob = pd.read_fwf(folder+'/'+n_id+'/'+'prob.txt', header = None)
        self.cost = pd.read_fwf(folder+'/'+n_id+'/'+'cost.txt', header = None)
        
        assert np.count_nonzero(self.cap) == len(self.prob), "Number of edges in the weighted adjacency matrix and the probability update matrix do not match."

        try:
            self.source, self.sink = pd.read_fwf(folder+'/'+n_id+'/'+'source_sink.txt', header = None).values[0]
        except: 
            #print('No source_sink.txt file found. Using the adjacency matrix to identify source and sink nodes.')
            self.source = self.cap.index[(self.cap.sum(axis = 0) == 0)].tolist()[0]
            self.sink = self.cap.index[(self.cap.sum(axis = 1) == 0)].tolist()[0]
            #print(f"Assumed nodes are \n Source node: {self.source} \n Sink node: {self.sink}")
        #generate the network 
        self.network = nx.from_pandas_adjacency(self.cap, create_using = nx.DiGraph())
        self.arc_order = list(self.network.edges())
        self.nEdges = len(self.network.edges())
        self.order_method = 'Default'
    def __str__(self):
        return f"Network id: {self.n_id} \nData is obtained from: {self.folder}"
    
    def orderArcsRandom(self, seed = 3):
        #set seed 
        np.random.seed(seed)
        #create indices to track shuffle locations
        indices = np.arange(0, self.nEdges)
        #shuffle 
        np.random.shuffle(indices)
        #reorder the arcs by the shuffled indices 
        self.arc_order = [x for _, x in sorted(zip(indices, self.arc_order))]
        #reorder rows of self.prob by the index 
        self.prob = self.prob.set_index(indices)
        self.order_method = 'Random'
         
        
    def orderArcsBetweennessCentrality(self, output = False):
        centrality = nx.edge_betweenness_centrality(self.network)
        def edge_betweenness_centrality(arc):
            return centrality[arc]
        def edge_betweenness_centrality_index(idx):
            return list(centrality.values())[idx]
        #create and sort the edge indices by centrality values 
        indices = list(np.arange(0, self.nEdges))
        indices.sort(key = lambda i: edge_betweenness_centrality_index[i])
        self.prob.index = indices
        #sort the arcs as well 
        self.arc_order.sort(key = lambda i: edge_betweenness_centrality[i])
        #record the sorting method 
        self.order_method = 'Betweenness Centrality'
        if output:
            return centrality 
        
    # def order_arcs_path_removal(self):
    #     arc_order = []   
    #     self.arc_order = arc_order 
    
    def create_scen_tree_v2(self, params):
        from timeit import default_timer as timer
        from treelib import Tree
        from math import prod as prod  
        import numpy as np
        from ortools.graph.python import max_flow
        
        # Set the parameters 
        G = self.network
        s = self.source 
        t = self.sink 
        R = params['budget']
        cutoff = params['cutoff']
        nEdges = self.nEdges
        self.mop_prob = self.prob.loc[:][0]

        class nodeData():
                def __init__(self):
                    self.capacities = None
                    self.solution_flows = None
                    self.mf = None
                    self.spent_budget = None 
                    self.scenario = None
                    self.mop = None
                    self.level = None
                    self.is_left = None
                    self.mf_lower = None 

        def check_resolve(arc, solution_flows):
            if solution_flows[arc] > 0:
                return True
            return False

        def remove_arc(arcs, capacities):
            new_capacities = capacities.copy()
            for arc in arcs:
                if self.mop_prob[arc] > 0:
                    new_capacities[arc] = 0 
            return new_capacities

        def mop_calc(scen, prob = self.mop_prob):
            scen_ub = prod([prob[arc] ** scen[1 + arc] for arc in range(len(scen[1:]))])
            return round(scen_ub, 4) 
              
        def branch_node(n_id):
            node_data = mytree.get_node(n_id).data
            #print(n_id)
            #print(f"mf: {node_data.mf}, r = {node_data.spent_budget}, mop: {node_data.mop}, level:{node_data.level}")
            #check if the node should be branched r pruned 
            if node_data.mf > 0 and node_data.spent_budget < R and node_data.mop > cutoff and node_data.level < nEdges:
                #Create the left children 
                left_id = len(mytree) + 1
                left_data = nodeData()
                left_data.capacities = remove_arc([node_data.level], node_data.capacities)
                left_data.level = node_data.level + 1
                #print(n_id, left_data.capacities)
                #check if resolve needed
                if check_resolve(node_data.level, node_data.solution_flows):
                    #remove the arc by setting the capacity of the corresponding arc to zero
                    #mf_time = timer()
                    smf = max_flow.SimpleMaxFlow()
                    all_arcs = smf.add_arcs_with_capacity(start_nodes, end_nodes, left_data.capacities)
                    status = smf.solve(s, t)
                    left_data.solution_flows = smf.flows(all_arcs)
                    left_data.mf = smf.optimal_flow()
                    #max_flow_time += timer() - mf_time
                    #mf_count += 1
                else:
                    left_data.solution_flows = node_data.solution_flows
                    left_data.mf = node_data.mf

                left_data.is_left = True
                left_data.spent_budget = node_data.spent_budget + 1 
                
                scen_list = node_data.scenario + [1]
                left_data.scenario = scen_list
                left_data.mop = mop_calc(scen_list)
                
                #tree_dict[left_id] = left_data
                mytree.create_node("n"+str(left_id),left_id, parent = n_id, data=left_data)
                branch_node(left_id)
                #create the right node 
                right_id = len(mytree) + 1 #node with no arc interdicted x_ij=0 
                right_data = nodeData()
                right_data.capacities = node_data.capacities
                right_data.solution_flows = node_data.solution_flows
                right_data.mf = node_data.mf
                right_data.spent_budget = node_data.spent_budget
                right_data.level = node_data.level + 1

                scen_list = node_data.scenario + [0]
                right_data.scenario = scen_list
                right_data.is_left = False
                right_data.mop = node_data.mop
                #tree_dict[right_id] = right_data
                mytree.create_node("n"+str(right_id),right_id, parent = n_id, data=right_data)
                branch_node(right_id)
            elif node_data.mop <= cutoff:
                 #remove the arc by setting the capacity of the corresponding arc to zero
                 #mf_time = timer()
                 smf = max_flow.SimpleMaxFlow()
                 caps = remove_arc(range(node_data.level, nEdges), node_data.capacities)
                 all_arcs = smf.add_arcs_with_capacity(start_nodes, end_nodes, caps)
                 status = smf.solve(s, t)
                 solution_flows = smf.flows(all_arcs)
                 node_data.mf_lower = smf.optimal_flow()                 
            else:        
                pass  
        
        #root node
        n_id = 1
        smf = max_flow.SimpleMaxFlow()
        #get start and end nodes 
        start_nodes, end_nodes = zip(*list(G.edges))
        start_nodes, end_nodes = list(start_nodes), list(end_nodes)
        capacities = [list(edge[2].values())[0] for edge in G.edges.data()]
        
        all_arcs = smf.add_arcs_with_capacity(start_nodes, end_nodes, capacities)
        smf.solve(s, t)
        solution_flows = smf.flows(all_arcs)
        mf = smf.optimal_flow()
                      
        #populate first node (use or-tools convention)
        nd = nodeData()
        nd.capacities = capacities
        nd.solution_flows = solution_flows
        nd.mf = mf
        nd.spent_budget = 0
        nd.scenario = [-1]
        nd.mop = 1 
        nd.is_left = 0
        nd.level = 0

        start = timer()        
        mytree = Tree()
        mytree.create_node("n"+ str(n_id),n_id,data=nd)
        branch_node(n_id)
        self.scen_tree_time = timer() - start  
        print('total seconds:', self.scen_tree_time)   
        print('# of leaves:', len(mytree.leaves()))                           
        
        self.scen_tree = mytree
        
    
    def create_scen_tree(self, params, rem_edge_list = None):
        from timeit import default_timer as timer
        from treelib import Tree
        import math 
        
        assert params['budget'] is not None, 'No budget key found in params dictionary'
        
        start = timer()
        # Create a tree object 
        mytree = Tree()
        G = self.network
        s = self.source 
        t = self.sink 
        
        #if no arc/layer removal order is given use the assigned order (default is edge list from networkx) 
        if rem_edge_list is None:
            rem_edge_list = self.arc_order

        n_id = 1
        #name, id, data(graph, max_flow, l/r (0,1) indicator, spent_budget, scenario, maximum prob that can get in)  
        mytree.create_node("n"+ str(n_id),n_id,data=(G,nx.maximum_flow(G,s,t,'weight')[0], 0, 0, [-1], 1))
        L = [1]
        p_id = 1

        R = params['budget']
    
        #branching into two nodes and checking the flows
        for rem_edge_idx in range(self.nEdges):
            print(f'Building the scenario tree current level: {rem_edge_idx + 1}/{self.nEdges}')
            if (len(L)>0):                   
                Ltemp = L.copy()
                for p_id in Ltemp: 
                    L.remove(p_id) #remove the node p_id as we are checking it now. 
                    #if max-flow is larger than 0 and budget spent upto now is less than budget continue checking
                    #and mytree.get_node(p_id).data[5] > params['cutoff']
                    node_data = mytree.get_node(p_id).data
                    if node_data[1]>0 and node_data[3]<R and round(node_data[5],5) > params['cutoff'] : #and budget cut
                        #create the left node
                        n_id += 1 #node with arc interdicted x_ij=1
                        parent_R = node_data[3] #cost spent until the parent
                        H = node_data[0].copy() #get the parent's graph 
                        rem_edge = rem_edge_list[rem_edge_idx] #remove the edge from the rem_edge_list
                        H.remove_edge(rem_edge[0],rem_edge[1]) 
                        scen_list = node_data[4].copy()
                        scen_list.append(1)
                        #add the node to the scenario tree with it's own values
                        
                        scen_probs = self.prob.loc[np.arange(0, len(scen_list[1:]))][0]
                        scen_ub = np.prod(np.power(scen_probs, scen_list[1:]))
                        
                        mytree.create_node("n"+str(n_id),n_id, parent = p_id, data=(H, nx.maximum_flow(H,s,t,'weight')[0],1, parent_R+1, scen_list, scen_ub))
                        L.append(n_id) #add it to the list to be checked 
                        
                        #create the right node 
                        n_id += 1 #node with no arc interdicted x_ij=0 
                        parent_G = node_data[0] #get parent's graph
                        scen_list = node_data[4].copy()
                        scen_list.append(0)
                        
                        scen_probs = self.prob.loc[np.arange(0, len(scen_list[1:]))][0]
                        scen_ub = np.prod(np.power(scen_probs, scen_list[1:]))
                        #remove the maximum flow with getting only the parent node - will halven the max flow calculations
                        mytree.create_node("n"+str(n_id),n_id, parent = p_id, data=(node_data[0], node_data[1],0, parent_R, scen_list, scen_ub))
                        L.append(n_id) #add it to the list to be checked 
                            
        self.scen_tree_time = timer() - start
        self.scen_tree = mytree
    
    #plot the original network with the source, sink and the 
    def plot(self):
        pos = nx.spring_layout(self.network)
        node_colors = dict(zip(list(self.network.nodes), ['white']*len(self.network.nodes)))
        node_colors[self.source] = 'red'
        node_colors[self.sink] = 'blue'
        nx.draw(self.network, pos, node_color = [node_colors[node] for node in self.network.nodes], with_labels = True, node_size = 500)

        #edge_labels = [for edge in self.network.edges]
        #nx.draw_networkx_edge_labels(self.network, pos, edge_labels = )

    def print_output(self, params):
        import json 
        import math
        total = sum([math.comb(self.nEdges, k) for k in range(params['budget'] + 1)])
        results = {'Network': self.n_id, 'Number of edges': self.nEdges,
                   'Scenario Clustering Threshold': params['cutoff'],
                   'Scenario Tree Generation Time (s)': round(self.scen_tree_time,2),
                   'Arc Order': self.order_method,
                   'Scenario Tree Number of Leaf Nodes': len(self.scen_tree.leaves()),
                   'Total Number of Scenarios under current budget': total, 
                   'Reduction Percentage (1 - # of Leaf/ # of Total Possible)': round(1 - len(self.scen_tree.leaves()) / total, 3),
                   'Total Number of Scenarios': 2 ** self.nEdges,
                   'Reduction Percentage (Total)': 1 - len(self.scen_tree.leaves()) / 2 ** self.nEdges}
        pretty_results = json.dumps(results, indent=4)
        print(pretty_results)

    def plotScenTree(self, params):
        assert self.scen_tree is not None, print('Create a scenario tree first')
        import graphviz
        self.scen_tree.to_graphviz("scen_tree.dot")
        import subprocess
        graph_name = self.n_id + "_" + str(params['budget']) + "_" +  str(params['cutoff']) + ".png"
        subprocess.call(["dot", "-Tpng", "scen_tree.dot", "-o", graph_name])


    def store_output(self, params):
        import pandas as pd
        import math 
        total = sum([math.comb(self.nEdges, k) for k in range(params['budget'] + 1)])
        df = pd.Series()
        df['network_id'] = self.n_id
        df['budget'] = params['budget']
        df['cutoff'] = params['cutoff']
        df['scenario tree time'] = round(self.scen_tree_time,2)
        df['number of scenarios(leaves)'] = len(self.scen_tree.leaves())
        df['scenario reduction within budget'] = round(1 - len(self.scen_tree.leaves()) / total, 3)
        df['scenario reduction total'] = 1 - len(self.scen_tree.leaves()) / 2 ** self.nEdges
        return df 
        
class Prob:
    import gurobipy as gp
    from gurobipy import GRB 
    
    def __init__(self, network, params):
        #all input parameters
        self.network = network

        #clustering paramters parameters 
        self.params = params 
   
        #gurobi parameters 
        self.grb_params = {'OutputFlag': 0, 'TimeLimit': 60 * 60 * 3}
   
    def __str__(self):
        return f"Problem generator for network: {self.network.n_id}"
    
    def arc_based_formulation(self, timelimit = 60*60*3):
        mytree = self.network.scen_tree 
        p = self.network.prob 
        budget = self.params['budget']
        depth = self.network.network.number_of_edges()
        T = self.network.prob.shape[1]
        
        len_N = len(mytree.nodes)
        N = mytree.nodes
        
        
        Time = range(T)
        
        arcs = range(mytree.depth())
        
        count_nodes = []
        edge_list = []
        for t in range(T-1):
            for k in range(t+1):
                count_nodes.append((t,k))
                edge_list.append([(t+1,k),(t+1,k+1)])
        adj_count = dict(zip(count_nodes, edge_list))
        G_count = nx.from_dict_of_lists(adj_count, create_using=nx.DiGraph)

        
        m = gp.Model()
        
        for par in self.grb_params:
            m.setParam(par, self.grb_params[par])
        
        i = m.addVars(len_N, T, ub = 1, name = 'in_flow')
        x = m.addVars(arcs, T, vtype = GRB.BINARY, name = 'interdiction')
        q = m.addVars(arcs, G_count.edges, name = 'Indicator') #indicator variable 

        #Budget 
        m.addConstrs(sum(x[arc,t] for arc in arcs) <= budget for t in Time)

        #Indicator Variable Constraints - to count the number of interdictions 
        for arc in arcs:
            #initial flow is 1 in all counting networks 
            m.addConstr(sum(q[arc,i,j] for (i,j) in G_count.out_edges((0,0))) == 1)
            #for each counting node 
            for cn in list(G_count.nodes()):
                #if there are outgoing arcs 
                if len(G_count.out_edges(cn)) >= 1:
                    #and if our counting node is not the initial node
                    if cn != (0,0):
                        #flow balance in = out 
                        m.addConstr(sum(q[arc,i,j] for (i,j) in G_count.in_edges(cn)) == sum(q[arc,i,j] for (i,j) in G_count.out_edges(cn)))
                    #i: no attempt arc(remain), j: attemp arc (increase) 
                    (ii,jj) = G_count.out_edges(cn)
                    m.addConstr(q[arc,ii[0],ii[1]] <= 1 - x[arc,ii[0][0]])
                    m.addConstr(q[arc,jj[0],jj[1]] <= x[arc,ii[0][0]])

        nodes = [n for n in mytree.nodes]
        leaf_nodes = [leave.identifier for leave in mytree.leaves()]
        nonleaf_nodes = list(set(nodes) - set(leaf_nodes))
        leaf_vals = [leave.data.mf for leave in mytree.leaves()]
        
        #initial flow is 1
        m.addConstrs(i[0,t] == 1 for t in Time)        
        
        def addNodeConstr(n,t):
            l_id = mytree.children(n)[0].identifier - 1
            r_id = mytree.children(n)[1].identifier - 1
            arc = mytree.depth(n)
            m.addConstrs(i[n-1, t] == i[l_id, t] + i[r_id, t] for t in Time)
            #no interdiction no split 
            
            m.addConstr(i[l_id, t] <= x[arc,t]) 
            if t==0:
                  m.addConstr(i[l_id,t] <= p[t][arc]*i[n - 1, t])
            elif t>=1: 
                for k in range(t+1):
                    m.addConstr(i[l_id,t] <= (p[k][arc]*i[n - 1, t] + 1 - sum(q[arc,i,j] for (i,j) in G_count.in_edges((t,k)))))

        for t in Time:
            [addNodeConstr(n,t) for n in nonleaf_nodes]        
        

        # #flow balance 
        # for n in nonleaf_nodes:
        #     l_id = mytree.children(n)[0].identifier - 1
        #     r_id = mytree.children(n)[1].identifier - 1
        #     arc = mytree.depth(n)
        #     m.addConstrs(i[n-1, t] == i[l_id, t] + i[r_id, t] for t in Time)
        #     #no interdiction no split 
            
        #     for t in Time:
        #         m.addConstr(i[l_id, t] <= x[arc,t]) 
        #         if t==0:
        #               m.addConstr(i[l_id,t] <= p[t][arc]*i[n - 1, t])
        #         elif t>=1: 
        #             for k in range(t+1):
        #                 m.addConstr(i[l_id,t] <= (p[k][arc]*i[n - 1, t] + 1 - sum(q[arc,i,j] for (i,j) in G_count.in_edges((t,k)))))
    
    
        m.setObjective(sum(leaf_vals[idx]*i[leaf_nodes[idx] - 1, t] for idx in range(len(leaf_nodes)) for t in Time))
        m.optimize()
             
        self.model = m
        self.formulation = 'Arc'
        self.interdictions = [x[s] for s in x if x[s].x > 0]
        self.q = q
        self.clean_solution_output()
    
        return m
        
    def path_based_formulation(self):
        import itertools
        from timeit import default_timer as timer

        #input 
        T = self.network.prob.shape[1]

        depth = self.network.network.number_of_edges()
        budget = self.params['budget']
        mytree = self.network.scen_tree
        p = self.network.prob
        
        #get leaf nodes as they are scenarios 
        leaves = mytree.leaves()
        leaves_ids = [leaf.identifier for leaf in leaves]
        path_vals = [leaf.data.mf for leaf in leaves]
        paths = mytree.paths_to_leaves()
        paths_ids = list(range(0, len(leaves)))


        #dict(node:[passing flows])
        from collections import defaultdict

        def create_node_path_dict(paths):
            node_path_dict = defaultdict(list)
        
            for path_id, path in enumerate(paths):
                for node in path:
                    node_path_dict[node].append(path_id)
        
            return dict(node_path_dict)
        start = timeit.default_timer()
        node_path_dict = create_node_path_dict(paths)
        end = timeit.default_timer()
        print('Paths are found, it took', round(end-start,2),'seconds.')        
                          
        Tk = list()
        for t in range(T):
            Tk.append(list(product([t], np.arange(t+1))))
        
        Tk = list(itertools.chain(*Tk)) 
        Tk.remove((0,0))

        Time = range(T)
        
        arcs = range(mytree.depth())
        
        count_nodes = []
        edge_list = []
        for t in range(T-1):
            for k in range(t+1):
                count_nodes.append((t,k))
                edge_list.append([(t+1,k),(t+1,k+1)])
        adj_count = dict(zip(count_nodes, edge_list))
        G_count = nx.from_dict_of_lists(adj_count, create_using=nx.DiGraph)

        # #upper bound for x - if prob is zero don't interdict
        x_UB = []
        for arc in arcs:
            if p[0][arc]==0:
                my_UB = 0
            else:
                my_UB = 1
            x_UB.append(repeat(my_UB,T))  

        start = timeit.default_timer()
        m = gp.Model() 
        for par in self.grb_params:
            m.setParam(par, self.grb_params[par])
              
        f = m.addVars(paths_ids, Time, name = 'Path')
        x = m.addVars(arcs, Time, vtype = GRB.BINARY, ub = x_UB, name = 'Interdiction') #interdiction
        q = m.addVars(arcs, G_count.edges, name = 'Indicator') #indicator variable 
        m.update()
        path_vars = [m.getVars()[0:len(f)][i::T] for i in range(T)]
        x_vars = [m.getVars()[len(f):len(f)+len(x)][i::T] for i in range(T)]

        #flow sum must be 1 
        for t in Time:
            m.addConstr(gp.LinExpr([1]*len(path_vars[t]), path_vars[t]) == 1)

        #Budget 
        m.addConstrs(sum(x[arc,t] for arc in arcs) <= budget for t in Time)

        #Indicator Variable Constraints
        for arc in arcs:
            #initial flow is 1 in all counting networks 
            m.addConstr(sum(q[arc,i,j] for (i,j) in G_count.out_edges((0,0))) == 1)
            #for each counting node 
            for cn in list(G_count.nodes()):
                #if there are outgoing arcs 
                if len(G_count.out_edges(cn)) >= 1:
                    #and if our counting node is not the initial node
                    if cn != (0,0):
                        #flow balance in = out 
                        m.addConstr(sum(q[arc,i,j] for (i,j) in G_count.in_edges(cn)) == sum(q[arc,i,j] for (i,j) in G_count.out_edges(cn)))
                    #i: no attempt arc(remain), j: attemp arc (increase) 
                    (i,j) = G_count.out_edges(cn)
                    m.addConstr(q[arc,i[0],i[1]] <= 1 - x[arc,i[0][0]])
                    m.addConstr(q[arc,j[0],j[1]] <= x[arc,i[0][0]])
                    
        def addNodeConstr(n, t):
            n = mytree[n]
            if n.data.is_left:
                node_id = n.identifier 
                p_id = mytree.parent(n.identifier).identifier
                node_depth = n.data.level
                #leftside make gurobi expression for node_id and parent_id
                node_path_vars = [path_vars[t][node] for node in node_path_dict[node_id]] 
                node_path_sum = gp.LinExpr([1]*len(node_path_vars), node_path_vars) 

                parent_path_vars = [path_vars[t][node] for node in node_path_dict[p_id]] 
                parent_path_sum = gp.LinExpr([1]*len(parent_path_vars), parent_path_vars) 

       
                m.addConstr(node_path_sum <= x[node_depth-1,t])
                if t==0:
                     m.addConstr(node_path_sum <= p[t][node_depth-1]*parent_path_sum)
                elif t>=1: 
                    for k in range(t+1):
                        m.addConstr(node_path_sum <= (p[k][node_depth-1]*parent_path_sum + 1 - sum(q[node_depth-1,i,j] for (i,j) in G_count.in_edges((t,k)))))
        
        start = timer()
        for t in Time:
             #in each period flow sum should be equal to 1 
             m.addConstr(gp.LinExpr([1]*len(path_vars[t]), path_vars[t]) == 1)
             [addNodeConstr(n,t) for n in mytree.nodes]
        print(f"time to build node constraints: {round(timer() - start, 2)}")
             
        m.setObjective(sum(gp.LinExpr(path_vals, path_vars[t]) for t in Time))
        end = timeit.default_timer()
        print('Starting optimization, it took', round(end - start, 2) ,'seconds to build the model')
        m.optimize()
        
        self.model = m
        self.formulation = 'Path'
        self.interdictions = [x[s] for s in x if x[s].x > 0]
        self.flows = [f[s] for s in f if f[s].x > 0] 
        self.q = q
        self.path_vars = path_vars 
        self.my_flow = node_path_dict
        self.clean_solution_output()

        return m
    
    def get_lower_bound(self):
        import itertools
        
        #TODO: if m exists modify the objective, else create the problem.
        
        #input 
        T = self.network.prob.shape[1]

        depth = self.network.network.number_of_edges()
        budget = self.params['budget']
        mytree = self.network.scen_tree
        p = self.network.prob
        
        #get leaf nodes as they are scenarios 
        leaves = mytree.leaves()
        leaves_ids = [leaf.identifier for leaf in leaves]
        path_vals = [leaf.data.mf_lower if leaf.data.mf_lower is not None else leaf.data.mf for leaf in mytree.leaves()]
        paths = mytree.paths_to_leaves()
        paths_ids = list(range(0, len(leaves)))


        #dict(node:[passing flows])
        from collections import defaultdict

        def create_node_path_dict(paths):
            node_path_dict = defaultdict(list)
        
            for path_id, path in enumerate(paths):
                for node in path:
                    node_path_dict[node].append(path_id)
        
            return dict(node_path_dict)
        start = timeit.default_timer()
        node_path_dict = create_node_path_dict(paths)
        end = timeit.default_timer()
        print('Paths are found, it took', round(end-start,2),'seconds.')        
                          
        Tk = list()
        for t in range(T):
            Tk.append(list(product([t], np.arange(t+1))))
        
        Tk = list(itertools.chain(*Tk)) 
        Tk.remove((0,0))

        Time = range(T)
        
        arcs = range(mytree.depth())
        
        count_nodes = []
        edge_list = []
        for t in range(T-1):
            for k in range(t+1):
                count_nodes.append((t,k))
                edge_list.append([(t+1,k),(t+1,k+1)])
        adj_count = dict(zip(count_nodes, edge_list))
        G_count = nx.from_dict_of_lists(adj_count, create_using=nx.DiGraph)

        # #upper bound for x - if prob is zero don't interdict
        x_UB = []
        for arc in arcs:
            if p[0][arc]==0:
                my_UB = 0
            else:
                my_UB = 1
            x_UB.append(repeat(my_UB,T))  

        start = timeit.default_timer()
        m = gp.Model() 
        for par in self.grb_params:
            m.setParam(par, self.grb_params[par])
              
        f = m.addVars(paths_ids, Time, name = 'Path')
        x = m.addVars(arcs, Time, vtype = GRB.BINARY, ub = x_UB, name = 'Interdiction') #interdiction
        q = m.addVars(arcs, G_count.edges, name = 'Indicator') #indicator variable 
        m.update()
        path_vars = [m.getVars()[0:len(f)][i::T] for i in range(T)]
        x_vars = [m.getVars()[len(f):len(f)+len(x)][i::T] for i in range(T)]

        #flow sum must be 1 
        for t in Time:
            m.addConstr(gp.LinExpr([1]*len(path_vars[t]), path_vars[t]) == 1)

        #Budget 
        m.addConstrs(sum(x[arc,t] for arc in arcs) <= budget for t in Time)

        #Indicator Variable Constraints
        for arc in arcs:
            #initial flow is 1 in all counting networks 
            m.addConstr(sum(q[arc,i,j] for (i,j) in G_count.out_edges((0,0))) == 1)
            #for each counting node 
            for cn in list(G_count.nodes()):
                #if there are outgoing arcs 
                if len(G_count.out_edges(cn)) >= 1:
                    #and if our counting node is not the initial node
                    if cn != (0,0):
                        #flow balance in = out 
                        m.addConstr(sum(q[arc,i,j] for (i,j) in G_count.in_edges(cn)) == sum(q[arc,i,j] for (i,j) in G_count.out_edges(cn)))
                    #i: no attempt arc(remain), j: attemp arc (increase) 
                    (i,j) = G_count.out_edges(cn)
                    m.addConstr(q[arc,i[0],i[1]] <= 1 - x[arc,i[0][0]])
                    m.addConstr(q[arc,j[0],j[1]] <= x[arc,i[0][0]])

        for t in Time:
            #in each period flow sum should be equal to 1 
            m.addConstr(gp.LinExpr([1]*len(path_vars[t]), path_vars[t]) == 1)
            for n in mytree.nodes:
                n = mytree[n]
                if n.data.is_left:
                    my_id = n.identifier 
                    p_id = mytree.parent(n.identifier).identifier #parent's id 
                    my_depth = mytree.depth(n)
                    m.addConstr(sum(f[i,t] for i in node_path_dict[my_id]) <= x[my_depth-1,t])
                    if t==0:
                         m.addConstr(sum(f[i,t] for i in node_path_dict[my_id]) <= p[t][my_depth-1]*sum(f[i,t] for i in node_path_dict[p_id]))
                    elif t>=1: 
                        for k in range(t+1):
                            m.addConstr(sum(f[i,t] for i in node_path_dict[my_id]) <= (p[k][my_depth-1]*sum(f[i,t] for i in node_path_dict[p_id]) + 1 - sum(q[my_depth-1,i,j] for (i,j) in G_count.in_edges((t,k)))))
    
        m.setObjective(sum(gp.LinExpr(path_vals, path_vars[t]) for t in Time))
        end = timeit.default_timer()
        print('Starting optimization, it took', round(end - start, 2) ,'seconds to build the model')
        m.optimize()
        
        self.lb = m.objVal
        
        return m
        
    def clean_solution_output(self):
        assert self.interdictions is not None, "Run a program before creating outputs"
        import re
        sols = []
        for sol in self.interdictions:
            input_string = sol.VarName
            # Define a regular expression pattern to match the value inside parentheses
            pattern = r"\[(\d+),(\d+)\]"
            # Use re.findall to find all occurrences of the pattern in the string
            matches = re.findall(pattern, input_string)
            values_inside_parentheses = matches[0]
            int_values = [int(value) for value in values_inside_parentheses]
            sols.append(int_values)
        
        d = dict()
        for s in sols:
            try: 
                d.get(s[1])
                d[s[1]].append(s[0])
            except:
                d[s[1]] = [s[0]]
        
        self.sols = d

        return d

    def store_output(self, params, df = pd.Series()):
        import pandas as pd 
        df['Formulation'] = self.formulation
        df['Found obj'] = round(self.model.objVal,5)
        df['Interdictions'] = self.sols    
        df['Optimization runtime'] = round(self.model.runtime,2)
        return df

#%% 
n_id = 'network_0'
n = Network(n_id)
params = {'budget':n.nEdges, 'cutoff': 0}
n.plot()
# n.create_scen_tree_v2(params)
# n.plotScenTree(params)



