from pyomo.opt import SolverStatus, TerminationCondition
from collections import Counter
from manage_graphs import compute_tree_size, extract_trees, build_graph, build_multitree_graph
import networkx as nx
from copy import deepcopy
from collections import defaultdict

donors = set()

ROUNDING_ERROR = 0.0005 # linearization includes some error

class WrongSolution(Exception):
    def __init__(self):
        super().__init__("The solution did not pass the tests:"
                         " Check the logs.")
    pass


def check_solution(model, results, g,  max_out_deg=0):
    error = False
    global donors
    donors = set() # this will be filled by check_l2_constraints
    print(donors)
    print('===========')
    print('Starting Tests')
    print('===========')
    gap = 0
    graph = build_graph(g, model, model.u, model.p)

    if (results.solver.status == SolverStatus.ok) and \
       (results.solver.termination_condition == TerminationCondition.optimal):
        print("Optimal Solution Found")
    elif (results.solver.status == SolverStatus.aborted):
        lb = results['Problem'][0]['lower bound']
        ub = results['Problem'][0]['upper bound']
        gap = (ub - lb)/ub
        print(ub,lb,gap)
        error = "Aborted"
    elif (results.solver.termination_condition == TerminationCondition.infeasible):
        print("ERROR: No feasible solution found!")
        error="Infeasible"
    else:
        # Something else is wrong
        print(f"ERROR: Solver Status: {results.solver.status}")
        error=True
    l2_err, max_diam = check_l2_constraints(model)
    flow_err = check_flow_constraints(model, results)
    if max_out_deg:
        tree_err = check_tree_size(model, graph, max_out_deg)
    else:
        tree_err = False
    error = error or l2_err or flow_err or tree_err
    
    #### TODO: test the corner case max_degree=0 and max_diam=0
    #model.f.pprint()
    if error:
        print('=======================')
        print('ERRORS are present')
        print('=======================')
    else:
        print('=======================')
        print('All tests passed')
        print('=======================')

    return len(donors), max_diam, error, gap, graph
    
def check_solution_multitree(model, results, g,  max_out_deg=0):
    error = False
    global donors
    donors = set()
    print('===========')
    print('Starting Tests')
    print('===========')
    gap = 0
    

    if (results.solver.status == SolverStatus.ok) and \
        (results.solver.termination_condition == TerminationCondition.optimal):
        print("Optimal Solution Found")
    elif (results.solver.status == SolverStatus.aborted):
        lb = results['Problem'][0]['lower bound']
        ub = results['Problem'][0]['upper bound']
        gap = (ub - lb)/ub
        print(ub,lb,gap)
        error = "Aborted"
    elif (results.solver.termination_condition == TerminationCondition.infeasible):
        print("ERROR: No feasible solution found!")
        error="Infeasible"
    else:
        # Something else is wrong
        print(f"ERROR: Solver Status: {results.solver.status}")
        error=True
    
    graph, tree_list = build_multitree_graph(g, model)
    tree_err = False
    multi_tree_donor_err = False
    l2_err = False
    max_diam = 0

    for tree in model.R:
        print(f"== Checking conditions on the multitree {tree} == ")
        l2_err_tree, diam = check_l2_constraints(model, tree)
        if diam > max_diam:
            diam = max_diam
        l2_err += l2_err_tree
        if max_out_deg:
            tree_err += check_tree_size(model, tree_list[tree-1], max_out_deg)

    donors_sets = defaultdict(set)
    for u in model.u:
        if model.u[u].value > ROUNDING_ERROR:
            if u[1] == 0:
                donors_sets[u[2]].add(u[0])

    for tree, d_set in donors_sets.items():
        if not d_set:
            print(f'ERROR: tree {tree} has no donors!')
            multi_tree_donor_err = True
        for right_tree, right_d_set in donors_sets.items():
            if right_tree == tree:
                continue
            if d_set & right_d_set:
                print(f'ERROR: tree {tree} and tree {right_tree} have shared',
                      f'donors: {d_set} and {right_d_set}')
                multi_tree_donor_err = True

    error = error or l2_err or multi_tree_donor_err or tree_err 
    if error:
        print('=======================')
        print('ERRORS are present')
        print('=======================')
    else:
        print('=======================')
        print('All tests passed')
        print('=======================')
    
    return len(donors), max_diam, error, gap, graph

def check_l2_constraints(model, tree_number=0, total_trees=1):
    """ tree_number will use only the nodes/edges of a specific tree in 
    the multitree configuration """
    distance = Counter()
    diam = max(model.D)
    if tree_number:
        model_u = deepcopy(model.u) 
        for u in list(model_u):
            if u[2] != tree_number:
                del model_u[u]
        model_p = deepcopy(model.p)
        for e in model_p:
            if e[2] != tree_number:
                del model_p[e]
    else:
        model_u = model.u
        model_p = model.p

    # check that diameter is respected
    max_diam = 0
    error = False
    for u in model_u:
        if model_u[u].value > ROUNDING_ERROR:
            if u[1] == 0:
                donors.add(u[0])
            if u[1] > max_diam:
                max_diam = u[1]
                if u[1] > diam:
                    print(f"ERROR: diameter of node {u[0]} is: {u[1]}> {diam}")
                    error = True
            distance[u[0]] += 1
            if distance[u[0]] > 1 and not tree_number:
                # if tree_number !=0 then we can have more than one distance
                print(f"ERROR: node {u[0]} has more than one distance variable")
                error = True
    inc_edges = Counter()
    out_edges = Counter()
    if not donors:
        print('ERROR: No Donors in the network!')
        error = True
    # check that incoming degree is respected
    for e in model_p:
        i = e[0] 
        j = e[1]
        
        if model_p[e].value > ROUNDING_ERROR: # some times solver gives a float 
                                       # solution very close to zero
            inc_edges[j] += 1
            out_edges[i] += 1
            if inc_edges[j] > total_trees:
                error = True
                print(f'ERROR: node {j} has {inc_edges[j]} incoming edges!')        
    
    # check only valid edges exist
    for e in model_p:
        i = e[0] 
        j = e[1]
        if model_p[e].value > ROUNDING_ERROR and not model.E[(i,j)]:
            print(f"ERROR: edge {(i,j)} was not present in the original graph")
            error=True
    # flow tests
    if not error:
        print(f'Diameter {diam} and single distance is respected')
        print('Incoming edges are respected')
        print('Only original edges are in the graph')
        print('At least one donor exists')

    return error, max_diam

def check_flow_constraints(model, results):
    error = False
    if model.find_component('f'):
        out_flow = {}
        in_flow = {}
        edge_flow = {}
        for i in model.V:
            out_flow[i] = 0  # donor nodes
            in_flow[i] = 0
        for i in model.V:
            for j in model.V:
                if (i,j) not in edge_flow:
                    edge_flow[(i,j)] = 0
                for h in model.V:
                    out_flow[i] += model.f[(i,j,h)].value                 
                    in_flow[j] += model.f[(i,j,h)].value   
                    edge_flow[(i,j)] += model.f[(i,j,h)].value
                    if model.p[(i,j)].value < ROUNDING_ERROR and \
                       model.f[(i,j,h)].value > ROUNDING_ERROR:
                        print(f"ERROR: edge {(i,j)} is off ",
                              f"and flow {(i,j,h)} is on")
                        error=True    
        total_demand = sum([d for (n,d) in model.demand.items() if\
                            n not in out_flow])
        tot_generated_flow = 0
        for i in out_flow:  
            if out_flow[i] > model.out_capacity[i] + ROUNDING_ERROR:
                print(f'ERROR: node {i} has more output flow than its ',
                      f'capacity {model.put_capacity[i]}')
                error=True
            if i in donors:
                tot_generated_flow += out_flow[i]
        if tot_generated_flow + ROUNDING_ERROR*len(donors)< total_demand:
            print(f'ERROR: total output from root nodes is  {tot_generated_flow} '
                  f'which is lower than the total demand {total_demand}')
            error=True
        for i in in_flow:
            if i in donors and in_flow[i] > ROUNDING_ERROR:
                print(f'ERROR: donor {i} has input flow {in_flow[i]}')
                error=True
            if i not in donors and in_flow[i] + ROUNDING_ERROR < model.demand[i]:
                print(f'ERROR: node {i} has flow {in_flow[i]}, less than demand {model.demand[i]}')
                error=True
        for (i,j),f in edge_flow.items():
            try:
                if f > model.lc[(i,j)] + ROUNDING_ERROR:
                    print(f'ERROR: link {(i,j)} has more flow than its capacity:',
                          f'{f} Vs {model.lc[(i,j)]}')
                    error=True
            except AttributeError:
                break
                
    if not error:
        print('Flow constraints are respected')
    return error
   
def check_tree_size(model, graph, max_out_deg):
    error = False
    tree_size = compute_tree_size(max(model.D.data()), max_out_deg)[max(model.D.data()), abs(max_out_deg)]
    filtered_g = extract_trees(graph)
    for c in nx.weakly_connected_components(filtered_g):
        if len(c) > tree_size:
            error = True
            print(f'ERROR: tree size is {len(c)} while it should be {tree_size}')
    if not error:
        print(f'Max tree size {tree_size} is respected')
    return error
        
    
    