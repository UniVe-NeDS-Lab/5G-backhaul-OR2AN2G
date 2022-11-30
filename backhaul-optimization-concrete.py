#! /usr/bin/env python3
import pyomo.environ as pyo
import sys
import os
import json
#from collections import Counter
import argparse 
from datetime import datetime
from tests import check_solution, check_solution_multitree, WrongSolution
from models import make_l2_model, make_flow_model, make_l2_multitree_model
import time
from manage_graphs import make_random_graph, save_graph, save_multitree_graph
from manage_graphs import make_poisson_graph, make_grid_graph, read_graph_file
from cell_score import  set_scenario, sample
import tables
from collections import defaultdict

        
def solve_model(model, args, solver='gurobi'):
    print('Solving the model')
    opt = pyo.SolverFactory(solver)
    if args.timelimit:
        if solver == 'cplex':
            opt.options['timelimit'] = int(args.timelimit*60)
        elif solver == 'glpk':         
            opt.options['tmlim'] = int(args.timelimit*60)
        elif solver == 'gurobi':           
            opt.options['TimeLimit'] = int(args.timelimit*60)   
    # maybe limit the running time
    results = opt.solve(model, logfile=run_path+"result.log") 
    # model.options['TimeLimit'] = 10
    with open(run_path+"model.log",'w') as f:
        out = sys.stdout
        sys.stdout = f
        model.pprint()
        sys.stdout = out
    return results
    

def parse_arguments():
    parser = argparse.ArgumentParser()
    gtype = parser.add_mutually_exclusive_group()
    gtype.add_argument('-n', type=int, help="Size of the graph")
    parser.add_argument('-r', type=int, help="Number or runs to repeat",
                        default=0)
    parser.add_argument('--type', choices=['erdosh', 'poisson', 
                                           'poisson_weighted', 'grid'], 
                        help="Use a synthetic graph of a certain kind", 
                        default='erdosh')
    gtype.add_argument('--graph', 
                        help="Provide a real graph in .graphml format")
    parser.add_argument('--l2model', choices=['single', 'multitree'], 
                        help="only one tree or multi-trees",
                        default='single')
    parser.add_argument('--flow', action='store_true',
                        help="consider flow allocation",
                        default=False)
    parser.add_argument('--diam', type=int, help="max diameter",default=3)
    parser.add_argument('--deg', type=int, default=0,
                        help="max outgoing degree (ignored with flow model)."
                        " If negative it starts with the absolute value and"
                        " decreases at each step")
    parser.add_argument('--capacity', help='maximum capacity of a node',
                        type=int, default=800)
    parser.add_argument('--seed', help='the random seed', type=int, default=0)
    parser.add_argument('--timelimit', help='limit the computation to X minutes',
                        type=float, default=0)

    args = parser.parse_args()
    if args.seed and args.r:
        parser.print_help()
        print("You can not set one seed for multuple runs: --seed and -r are"
              " incompatible")
        exit()
    if not args.n and not args.graph:
        parser.print_help()
        print('Provide at least one between --graph or -n')
        exit()
    return args
    

def set_up_run(args):

    demand=100
    capacity = args.capacity
    seed = -1
    if args.graph:
        g = read_graph_file(args.graph)
    elif args.type == 'grid':
        g, seed = make_grid_graph(args.n)
    elif args.type == 'erdosh':
        g, seed = make_random_graph(args.n, seed=args.seed)
    elif args.type == 'poisson':
        g, coords, seed = make_poisson_graph(args.n, seed=args.seed)
    elif args.type == 'poisson_weighted':
        g, coords, seed = make_poisson_graph(args.n, seed=args.seed)
        circles = set_scenario(radius=100, coords=coords)
        circle_score, _ = sample(circles)
        max_score = 0
        for i in circle_score:
            g.nodes[i]['weight'] = circle_score[i] 
            #g.nodes[i]['x'] = coords[i][0]
            #g.nodes[i]['y'] = coords[i][1]
            max_score = max(max_score, circle_score[i])
        demand=circle_score
        if max_score > args.capacity:
            print(f'ERROR: Capacity ({args.capacity}) is lower than maximum'
                  f' weight ({max_score}), please use a resonable number')
            exit()
    else:
        print('ERROR, unknown graph type {graph_type} and no file provided')
        exit()

    model = pyo.ConcreteModel()
    if not args.deg:
        print("The unerlying l2 model needs a degree parameter")
        exit()
    if args.l2model != "multitree":
        make_l2_model(model, nodes=list(g.nodes()), edges=list(g.edges()), 
                      max_out_deg=args.deg, max_diam=args.diam)
    else:
        make_l2_multitree_model(model, nodes=list(g.nodes()),
                                edges=list(g.edges()), 
                                max_out_deg=args.deg, max_diam=args.diam)
    if args.flow:
        link_capacity=defaultdict(dict)
        mcs_table = tables.MCS_distance_table()
        for edge in g.edges(data=True):
            if "distance" in edge[2]:
                bitrate = mcs_table.get_bitrate(edge[2]['distance'])
                if bitrate:
                    link_capacity[edge[0]][edge[1]] = bitrate
        #FIXME what is the maximum outgoing capacity per radio?
        make_flow_model(model, capacity_list=capacity, demand_list=demand, 
                        link_capacity=link_capacity)
    r = solve_model(model, args)
    return model, r, g, seed




args = None
result_path = 'results/'
run_path = ''
os.makedirs(result_path, exist_ok=True)

if __name__ == "__main__":
    args = parse_arguments()
    if args.r:
        results = {}
        result_path += f'{args.type}_{args.l2model}_{args.n}-'
        timestring = datetime.now().strftime("%d-%m-%H:%M:%S/")
        result_path += timestring
        os.makedirs(result_path)
        results['config'] = vars(args)
        with open(result_path+'config.json', "w") as f:
            json.dump(results,f) 
        results['runs'] = []
        for run in range(args.r):
            start_time = time.time()
            run_path = result_path+f'{run}/'
            os.makedirs(run_path)
            print(f'===================== Run ID: {run}')
            model, r, g, seed = set_up_run(args)
            if args.l2model != "multitree":
                donors, max_diam, error, gap, newg = check_solution(model, r, g, 
                                                          max_out_deg=args.deg)
            else:
                donors, max_diam, error, gap, newg = check_solution_multitree(model, 
                                                          r, g, 
                                                          max_out_deg=args.deg)                
            res = {'run_id':run, 'donors':donors, 'diameter':max_diam,
                   'time':f'{int(time.time() - start_time)}s', 
                   'error_status':error, 'seed':seed, 'gap':gap}
            results['runs'].append(res)
            if error and not gap:
                with open(result_path+'results.json', "w") as f:
                    json.dump(results,f)
                    save_graph(newg, outfolder=run_path)
                raise WrongSolution
            if args.l2model != "multitree":
                save_graph(newg, outfolder=run_path)
            else:
                save_multitree_graph(newg, outfolder=run_path)
       
        with open(result_path+'results.json', "w") as f:
                json.dump(results,f)

    else: # just one run
        run_path = result_path
        model, r, g, _= set_up_run(args)
        if args.l2model != "multitree":
            donors, max_diam, error, gap, newg = check_solution(model, r, g, 
                                                      max_out_deg=args.deg)
        else:
            donors, max_diam, error, gap, newg = check_solution_multitree(model, 
                                                      r, g, 
                                                      max_out_deg=args.deg)
        if args.l2model != "multitree":
            save_graph(newg, outfolder=run_path, labels=False)
        else:
            save_multitree_graph(newg, outfolder=run_path, labels=False)        
        if error and not gap:
            raise WrongSolution
 
