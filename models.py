#!/usr/bin/env python3

import pyomo.environ as pyo


def make_l2_multitree_model(model, nodes, edges, max_diam, max_out_deg, 
                            num_trees=2):
    
    model.D = pyo.RangeSet(0, max_diam)
    model.R = pyo.RangeSet(1, num_trees)
    model.D_pos = pyo.RangeSet(1, max_diam)
    model.V = pyo.Set(initialize=nodes)

    def edge_init(model, i, j):
        return (i,j) in edges
    
    model.E = pyo.Param(model.V, model.V, initialize=edge_init, 
                        domain=pyo.Binary)
    model.u = pyo.Var(model.V, model.D, model.R, domain=pyo.Binary, initialize=0)
    model.p = pyo.Var(model.V, model.V, model.R, domain=pyo.Binary, initialize=0)
    
    def objective(model):
        return sum(model.u[i,0,k] for i in model.V for k in model.R)
    
    def dist(model,i,k):
        return sum([model.u[i,l,k] for l in model.D]) == 1 - \
               sum(model.u[i,0,r] for r in model.R if k!=r)

    def single_root(model,i):
        return sum(model.u[i,0,k] for k in model.R) <= 1
    
    def one_incoming(model, j, k):
        return sum([model.p[i,j,k] for i in model.V]) == 1 - \
               sum(model.u[j,0,r] for r in model.R)
    
    def max_deg(model, i, k):
        return sum([model.p[i,j,k] for j in model.V]) <= model.deg
        
    def max_deg_decr(model, i, l, k):
        return sum([model.p[i,n] for n in model.V]) <= model.deg - l*model.u[i,l] 
 
    def incremental_distance(model, i, j, l, k):
        return model.p[i,j,k] <= 1 - model.u[j,l,k] + model.u[i,l-1,k]
    
    def only_existing_edges(model, i, j, k):
        return model.p[i,j,k] <= model.E[i,j]
  
    def only_one_direction(model, i, j, k):
        return model.p[i,j,k] + model.p[j,i,k] <= 1
    
    def separate_trees(model, i, j):
        return sum(model.p[i,j,k] + model.p[j,i,k] for k in model.R) <= 1
    
    print('Building the L2 Multigraph model')
    
    model.OBJ = pyo.Objective(rule=objective)
    
    model.distance = pyo.Constraint(model.V, model.R, rule=dist)
    model.single_root = pyo.Constraint(model.V, rule=single_root)
    model.one_incoming = pyo.Constraint(model.V, model.R, rule=one_incoming)
    if max_out_deg > 0:
        model.deg = max_out_deg
        model.max_deg = pyo.Constraint(model.V, model.R, rule=max_deg)
    elif max_out_deg < 0: 
        model.deg = -max_out_deg
        model.max_deg = pyo.Constraint(model.V, model.D, model.R, 
                                       rule=max_deg_decr)

    model.incremental_distance = pyo.Constraint(model.V, model.V, model.D_pos, 
                                       model.R, rule=incremental_distance)
    model.only_existing_edges = pyo.Constraint(model.V, model.V, model.R,
                                       rule=only_existing_edges)
    model.unidirectional = pyo.Constraint(model.V, model.V, model.R,
                                       rule=only_one_direction)    
    model.separate_trees = pyo.Constraint(model.V, model.V, 
                                       rule=separate_trees) 
    
    return model


def make_l2_model(model, nodes, edges, max_diam, max_out_deg=0):
    model.D = pyo.RangeSet(0, max_diam)
    model.D_pos = pyo.RangeSet(1, max_diam)
    model.V = pyo.Set(initialize=nodes)
    
    def edge_init(model, i, j):
        return (i,j) in edges
    
    model.E = pyo.Param(model.V, model.V, initialize=edge_init, 
                        domain=pyo.Binary)
    model.u = pyo.Var(model.V, model.D, domain=pyo.Binary, initialize=0)
    model.p = pyo.Var(model.V, model.V, domain=pyo.Binary, initialize=0)
 
    
    def objective(model):
        return sum([model.u[i,0] for i in model.V])
    
    def dist_c(model,i):
        return sum([model.u[i,l] for l in model.D]) == 1
    
    def one_incoming(model, j):
        return sum([model.p[i,j] for i in model.V]) == 1 - model.u[j,0]
    

    def max_deg(model, i):
        return sum([model.p[i,j] for j in model.V]) <= model.deg
    
    def max_deg_decr(model, i, j):
        return sum([model.p[i,k] for k in model.V]) <= model.deg - j*model.u[i,j] 
    
    def incremental_distance(model, i, j, l):
        return model.p[i,j] <= 1 - model.u[j,l] + model.u[i,l-1]
        
    def distance_implies_incoming(model, j):
        return sum([model.u[j,l] for l in model.D_pos]) <= sum([model.p[i,j] for i in model.V])
    
    def only_existing_edges(model, i, j):
        return model.p[i,j] <= model.E[i,j]
    
    def only_one_direction(model, i, j):
        return model.p[i,j] + model.p[j,i] <= 1
   
    print('Building the L2 model')
    
    model.OBJ = pyo.Objective(rule=objective)
    
    model.distance = pyo.Constraint(model.V, rule=dist_c)
    model.one_incoming = pyo.Constraint(model.V, rule=one_incoming)
    if max_out_deg > 0:
        model.deg = max_out_deg
        model.max_deg = pyo.Constraint(model.V, rule=max_deg)
    elif max_out_deg < 0: 
        model.deg = -max_out_deg
        model.max_deg = pyo.Constraint(model.V, model.D, rule=max_deg_decr)

    model.incremental_distance = pyo.Constraint(model.V, model.V, model.D_pos, 
                                       rule=incremental_distance)
    model.distance_incoming = pyo.Constraint(model.V, rule=distance_implies_incoming)
    model.only_existing_edges = pyo.Constraint(model.V, model.V, 
                                       rule=only_existing_edges)
    model.unidirectional = pyo.Constraint(model.V, model.V, 
                                       rule=only_one_direction)    
    return model


def make_flow_model(model, capacity_list, demand_list, 
                    link_capacity={}):
 
    # if capacity_dict and demand_list are lists with 
    # nodes as indexing keys, this still works
    model.out_capacity = pyo.Param(model.V, initialize=capacity_list, 
                               domain=pyo.PositiveIntegers)
    if link_capacity:
        def fill_link_cap(model,i,j):
            try:
                cap = link_capacity[i][j]
            except KeyError:
                cap = 0
            return cap
        model.lc = pyo.Param(model.V, model.V,
                             initialize=fill_link_cap, 
                             domain=pyo.NonNegativeIntegers)
    model.demand = pyo.Param(model.V, initialize=demand_list,
                             domain=pyo.NonNegativeIntegers)
    model.f = pyo.Var(model.V, model.V, model.V, 
                      domain=pyo.NonNegativeIntegers, initialize=0)
    
    
    def capacity_limit(model, i):
        return sum([model.f[i,j,h] for j in model.V for h in model.V]) \
               <= model.out_capacity[i]
    
    def flow_conservation(model, i):
        return sum([model.f[i,j,h] for j in model.V for h in model.V if i!=h])\
               - \
               sum([model.f[j,i,h] for j in model.V for h in model.V if i!=h])\
               == model.u[i,0]*model.out_capacity[i]
               
    def flow_demand(model, j):
        return sum([model.f[i,j,j] for i in model.V]) >= \
               (1-model.u[j,0])*model.demand[j]
    
    def no_outgoing_flow(model, i, j):
        return model.f[i,j,i] == 0
    
    def selected_edges_only(model, i, j, h):
        return model.f[i,j,h] <= model.p[i,j]*model.out_capacity[i]
    
    def respect_link_capacity(model, i, j):
        return sum(model.f[i,j,h] for h in model.V) <= model.lc[i,j]

    print('Building the flow model')

    model.capacity_limit = pyo.Constraint(model.V, rule=capacity_limit)
    model.flow_conservation = pyo.Constraint(model.V, rule=flow_conservation)
    model.flow_demand = pyo.Constraint(model.V, rule=flow_demand)
    model.no_outgoing = pyo.Constraint(model.V, model.V, 
                                        rule=no_outgoing_flow)
    model.only_selected_edges = pyo.Constraint(model.V, model.V, model.V, 
                                        rule=selected_edges_only)
    if link_capacity:
        model.respect_link_capacity = pyo.Constraint(model.V, model.V, 
                                                     rule=respect_link_capacity)
    return model
