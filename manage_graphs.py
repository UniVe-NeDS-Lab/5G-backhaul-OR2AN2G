import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import lognorm
from numpy import exp, log, pi
from numpy.random import uniform
from numpy.random import seed as np_seed
import random 
import datetime

def read_graph_file(f):
    g = nx.read_graphml(f)
    for e in list(g.edges(data=True)):
        if e[2]['type'] == 'wired':
            g = nx.contracted_nodes(g, e[0], e[1], self_loops=False)
            del g.nodes[e[0]]['contraction']
    for e in list(g.edges(data=True)):
        if e[2]['los'] != True:
            g.remove_edge(e[0], e[1])
    
    return g.to_directed()

def make_line(n):
    g = nx.path_graph(n).to_directed()
    return g

def set_seed(seed):
    if not seed:
        seed = int(datetime.datetime.now().timestamp())
    random.seed(seed)
    np_seed(seed)
    print(f'using seed {seed}')
    return seed
        
def make_grid_graph(n):
    g = nx.grid_2d_graph(n,n).to_directed()
    newg = nx.DiGraph()
    for n in g:
        print(n)
        newg.add_node(f'{n[0]}-{n[1]}', x=n[0], y=n[1])
    for e in g.edges():
        newg.add_edge(f'{e[0][0]}-{e[0][1]}', f'{e[1][0]}-{e[1][1]}' )
    print(g.nodes())
    print(newg.nodes())
    return newg, 0

def make_random_graph(n, prob=0, seed=0):
    seed = set_seed(seed)
    if not prob:
        prob = log(n)/n
    for i in range(100):
        g = nx.erdos_renyi_graph(n, prob).to_directed()
        if nx.is_strongly_connected(g):
            break
    else:
        print("Disconnected graph, increase the edge probability")
        exit()
    print(f'Created a graph with {len(g)} nodes and {len(g.edges())} directed archs')
    return g, seed

def make_poisson_graph(n, density=50, seed=0, border=100):
    # density = nodes per squared km
    area = n/density
    size = int(area**0.5*1000) + border*2
    coords = {}
    sigma=0.9
    mu=3.04,
    seed = set_seed(seed)
    F = lognorm(s=sigma, scale=exp(mu)).cdf
    def plos(F, d, h=9, l=1/1430, r=8.07):
       return exp(-l*(1-F(h))*(2*r*d-pi*r**2))
   
    found = False
    while not found:
        for i in range(100):
            g = nx.DiGraph()
            for i in range(1,n+1):
                coords[i] = [uniform(border,size-border),
                             uniform(border,size-border)]
                g.add_node(i, x=coords[i][0], y=coords[i][1])
            for i in range(1,n+1):
                for j in range(i+1,n+1):
                    d = ((coords[i][0]-coords[j][0])**2 + 
                         (coords[i][1]-coords[j][1])**2)**0.5
                    if uniform() < plos(F,d):
                        g.add_edge(i,j)
                        g.add_edge(j,i)
            if nx.is_strongly_connected(g):
                found = True
                break
        else:
            density = density*1.1
            area = n/density
            size = int(area**0.5*1000) + border*2
            print(f"Can not find a connected graph, increasing the density to {density} ")
            
    print(f'Created a graph with {len(g)} nodes and {len(g.edges())}'
          f' directed archs on a {size}x{size}m area, with density {density}')
    return g, coords, seed

def compute_tree_size(max_diam, max_out_deg, decr=0):
    tree_size = {}
    if max_out_deg < 0:
        max_out_deg = -max_out_deg
        decr = 1
    for j in range(2, max_out_deg+1):
        for i in range(1, max_diam+1):
            nodes = 1
            level_nodes = 1
            degree = j
            for r in range(1, i+1):
                if not degree:
                    break
                level_nodes = level_nodes*degree
                degree -= decr
                nodes += level_nodes
            tree_size[(i,j)] = nodes
    return tree_size


def build_graph(g, model, model_u, model_p):
    donors_count = 0
    newg = nx.DiGraph()
    newg.add_nodes_from(g.nodes(data=True))
    for u in model_u:
        if model_u[u].value >= 0.5:
            if u[1] == 0:
                # then donor
                newg.nodes[u[0]]['label'] ='donor'
                donors_count += 1
            else:
                newg.nodes[u[0]]['label'] ='gNB'
    for e in model_p:
            if model_p[e].value > 0.5:
                f = 0
                if model.find_component('f'):
                    f = sum([model.f[e[0], e[1], h].value for h in model.V]) 
                newg.add_edge(e[0], e[1], label=f'{f}')
            elif e in g.edges():
                newg.add_edge(e[0], e[1], label='off')
    newg.nodes()
    return newg

    
def build_multitree_graph(g, model):
    multitree = []
    newg = nx.DiGraph()
    newg.add_nodes_from(g.nodes(data=True))
    for tree in model.R:
        sub_tree = nx.DiGraph()
        for u in model.u:
            if model.u[u].value >= 0.5:
                if u[1] == 0:
                    # then donor
                    newg.nodes[u[0]]['label'] =f'donor_{u[2]}'
                    if u[2] == tree:
                        sub_tree.add_node(u[0], label=f'donor_{u[2]}')
                else:
                    newg.nodes[u[0]]['label'] ='gNB'
                    if u[2] == tree:
                        sub_tree.add_node(u[0] ,label='gNB')
                
        for e in model.p:
            if e[2] != tree:
                continue
            if (e[0], e[1]) in g.edges():  # only existing edges
                if (e[0], e[1]) not in newg.edges():
                    if model.p[e].value > 0.5:
                        newg.add_edge(e[0], e[1], label=f'{e[2]}') 
                        sub_tree.add_edge(e[0], e[1], label='on')
                    else:
                        newg.add_edge(e[0], e[1], label='off') 
                else:
                    if model.p[e].value > 0.5:
                        newg[e[0]][e[1]]['label']=f'{e[2]}'
        multitree.append(sub_tree)
    return newg, multitree
                        
def extract_trees(g):
    def is_on(fr, to, *args):
        return g[fr][to]['label'] not in ['off']
    return nx.subgraph_view(g, filter_edge=is_on)    



def set_node_colors(newg, donors_main, gNBs):
    pos = {}
    try:
        for n in newg.nodes(data=True):
            pos[n[0]] = (n[1]['x'], n[1]['y'])
    except KeyError:
        try:
            pos = nx.nx_agraph.graphviz_layout(newg)
        except ImportError:
            pos = nx.kamada_kawai_layout(newg)
    try:
        donors_color = []
        gnb_color = []
        max_weight = max([w[1]['weight'] for w in newg.nodes(data=True)])
        min_weight = min([w[1]['weight'] for w in newg.nodes(data=True)])

        for n in donors_main:
            donors_color.append(newg.nodes[n]['weight'])    
        for n in gNBs:
            gnb_color.append(newg.nodes[n]['weight'])

        color_map = matplotlib.cm.ScalarMappable(cmap="Reds", 
                    norm=matplotlib.colors.Normalize(vmin=min_weight, 
                                                     vmax=max_weight,))

    except KeyError:
        gnb_color = '#0000ff'
        donors_color = '#ff0000'
        color_map = {}
    return pos, gnb_color, donors_color, color_map

def save_multitree_graph(newg, edges_label=True, outfolder='./', labels=True):
    
    nx.write_graphml_xml(newg, outfolder+'graph.graphml')
    f = plt.figure()

    all_donors = [n[0] for n in newg.nodes(data=True) if 'donor' in n[1]['label']]     
    gNBs = [n[0] for n in newg.nodes(data=True) if n[1]['label'] == 'gNB'] 
    
    edges_on = {(e[0], e[1]):int(e[2]['label']) for e in newg.edges(data=True) if 
                     e[2]['label'] != 'off'}
    edges_off = [(e[0], e[1]) for e in newg.edges(data=True) 
                     if e[2]['label'] == 'off'] 
    
    print('donors:', all_donors)
    print('gNB:', gNBs)
    pos,gnb_color,_,_ = set_node_colors(newg, all_donors, gNBs)
    # we use only the positions and gnb color
    node_size=300

    if labels:
        nx.draw_networkx_labels(newg, pos)
    else:
        node_size=100
        
    donor_color_list = []
    for d in all_donors:
        color = 1   # there can be donors without children, we set them to the
                    # color of tree 1
        for e in newg.out_edges(d, data=True): # is this a donor with children?
            if e[2]['label'] != 'off':         # find the tree it belongs to
                color = int(e[2]['label'])
                break
        donor_color_list.append(color)
            
    # draw gNBs
    nx.draw_networkx_nodes(newg, pos, nodelist=gNBs, node_shape='o',
                           node_color=gnb_color, node_size=node_size)
    # draw donors with colors
    nx.draw_networkx_nodes(newg, pos, nodelist=all_donors, node_shape='s',
                           node_color=donor_color_list, 
                           cmap=plt.cm.winter, linewidths=1, 
                           edgecolors="Black", 
                           node_size=node_size, vmin=1, 
                           vmax=max(donor_color_list))
  
    # draw off edges
    nx.draw_networkx_edges(newg, pos, edgelist=edges_off, style=':',
                           arrows=False, alpha=0.1) 

    
    edge_color_list = list(edges_on.values())

    color_vmax = 1 if not edge_color_list else max(edge_color_list)
    # draw on edges

    nx.draw_networkx_edges(newg, pos, edgelist=list(edges_on.keys()), 
                           connectionstyle='arc3, rad = 0.1', 
                           edge_color=edge_color_list, edge_cmap=plt.cm.winter,
                           edge_vmin=1, edge_vmax=color_vmax)
  
    plt.title(f'Nodes: {len(newg)}; Archs (directed):{len(newg.edges())};'
              f' Active Archs: {len(edges_on)}; Donors: {len(all_donors)}')
    f.savefig(outfolder+'graph.png', dpi=900)
    plt.close()

def save_graph(newg, edges_label=True, outfolder='./', labels=True):
    
    nx.write_graphml_xml(newg, outfolder+'graph.graphml')
    f = plt.figure()
    donors = [n[0] for n in newg.nodes(data=True) if n[1]['label'] == 'donor'] 
    gNBs = [n[0] for n in newg.nodes(data=True) if n[1]['label'] == 'gNB'] 
    edges_on = {(e[0], e[1]):e[2]['label'] for e in newg.edges(data=True) if 
                e[2]['label'] != 'off'}
    edges_labels = {(e[0], e[1]):e[2]['label'] for e in newg.edges(data=True) if 
                e[2]['label'] not in ['off', '0']}
    edges_off = [(e[0], e[1]) for e in newg.edges(data=True) if e[2]['label'] == 'off'] 
    print('donors:', donors)
    print('gNB:', gNBs)
    pos, gnb_color, donors_color, color_map = set_node_colors(newg, donors, gNBs)
  
    node_size=300
    if labels:
        nx.draw_networkx_labels(newg, pos)
    else:
        node_size=100
        
    nx.draw_networkx_nodes(newg, pos, nodelist=gNBs, node_shape='o',
                           node_color=gnb_color, cmap='Reds', node_size=node_size)
    nx.draw_networkx_nodes(newg, pos, nodelist=donors, node_shape='s',
                           node_color=donors_color, cmap='Reds', 
                           linewidths=1, edgecolors="Black", node_size=node_size)
    if color_map:
        plt.colorbar(color_map)
    nx.draw_networkx_edges(newg, pos, edgelist=edges_off, style=':',
                           arrows=False, alpha=0.1) 
    nx.draw_networkx_edges(newg, pos, edgelist=edges_on.keys())
    
    if edges_label:
        nx.draw_networkx_edge_labels(newg, pos, edge_labels=edges_labels,
                                     font_size=7)                      
    plt.title(f'Nodes: {len(newg)}; Archs (directed):{len(newg.edges())};'
              f' Active Archs: {len(edges_on)}; Donors: {len(donors)}')
    f.savefig(outfolder+'graph.png', dpi=900)
    plt.close()
    