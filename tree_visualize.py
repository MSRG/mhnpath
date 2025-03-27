import json
import matplotlib.pyplot as plt
import networkx as nx
from rdkit import Chem
from rdkit.Chem import Draw
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from reaction_cond import get_solvent_reagent

class Node:
    def __init__(self, smiles, cost_usd_per_g, depth):
        self.smiles = smiles
        self.cost_usd_per_g = cost_usd_per_g
        self.depth = depth
        self.subtrees = []

class Edge:
    def __init__(self, reaction_smiles, temperature, score, enzyme, rule, label):
        self.reaction_smiles = reaction_smiles
        self.temperature = temperature
        self.enzyme = enzyme
        self.score = score
        self.rule = rule
        self.label = label
        self.solvent = None
        self.reagent = None

def json_to_tree(json_file):
    try:
        with open(json_file, 'r') as file:
            data = json.load(file)
    except FileNotFoundError:
        print(f"Error: File '{json_file}' not found.")
        return None
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON file '{json_file}': {e}")
        return None
    
    if not isinstance(data, dict):
        print(f"Error: JSON data in '{json_file}' should be a dictionary.")
        return None
    
    try:
        root_node = build_node(data, True)
    except KeyError as e:
        print(f"Error: Missing key '{e}' in the JSON data.")
        return None
    
    return root_node

def build_node(node_data, root):
    smiles = node_data['smiles']
    cost_usd_per_g = node_data['cost_usd_per_g']
    depth = node_data['depth']
    
    node = Node(smiles, cost_usd_per_g, depth)
    
    if 'subtrees' in node_data and isinstance(node_data['subtrees'], list):
        for subtree_data in node_data['subtrees']:
            edge = build_edge(subtree_data)
            subtree_node = build_node(subtree_data['subtree'], True)
            if root:
                if node.subtrees != []:
                    rxn = node.subtrees[0][0].reaction_smiles
                    if edge.reaction_smiles == rxn:
                        node.subtrees.append((edge, subtree_node))
                else:
                    node.subtrees.append((edge, subtree_node))
            else:
                node.subtrees.append((edge, subtree_node))
    
    return node

def build_edge(edge_data):
    reaction_smiles = edge_data['reaction_smiles']
    temperature = edge_data['temperature']
    enzyme = edge_data['enzyme']
    score = edge_data['score']
    rule = edge_data['rule']
    label = edge_data['label']
    type_dis = edge_data['type_dis']
    try:
        solvent, reagent = get_solvent_reagent(reaction_smiles)
    except:
        solvent = None
        reagent = None
    
    edge = Edge(reaction_smiles, temperature, enzyme, score, rule, label, type_dis)
    edge.solvent = solvent
    edge.reagent = reagent
    return edge

# Function to draw a node with the 2D image of the molecule and cost
def draw_node(ax, node, pos, index):
    print(node.smiles)
    m = Chem.MolFromSmiles(node.smiles)
    img = Draw.MolToImage(m, size=(220, 220))
    imgbox = OffsetImage(img, zoom=0.4)
    imgbox.image.axes = ax
    
    ab = AnnotationBbox(imgbox, pos[index], frameon=False)
    ax.add_artist(ab)
    
    cost_text = f"Cost: ${round(node.cost_usd_per_g, 2)}/g"
    ax.annotate(cost_text, xy=pos[index], xytext=(0, -40), textcoords='offset points', ha='center', fontsize=8, color='red')

# Function to draw an edge with temperature, score, and type_dis only if it's from the root node
def draw_edge(ax, edge, start_pos, end_pos, color, is_from_root):
    mid_pos = ((start_pos[0] + end_pos[0]) / 2, (start_pos[1] + end_pos[1]) / 2)
    if is_from_root:
        edge_text = f"T={round(edge.temperature, 2)}°C, Dis={edge.type_dis}\n Sol={edge.solvent}, Reag={edge.reagent}"
    else:
        edge_text = f"T={round(edge.temperature, 2)}°C\n Sol={edge.solvent}, Reag={edge.reagent}"
    ax.annotate(edge_text, xy=mid_pos, xytext=(0, 10), textcoords='offset points', ha='center', fontsize=8, color=color)

def visualize_tree(root_node):
    G = nx.DiGraph()
    pos = {}
    labels = {}
    edge_colors = {}
    edge_color_mapping = {}
    current_index = 0

    # Traverse the tree and add nodes and edges to the graph
    def traverse(node, parent_index=None, parent_edge=None, is_from_root=False):
        nonlocal current_index
        node_index = current_index
        current_index += 1

        pos[node_index] = (node.depth, -node_index)
        labels[node_index] = node.smiles
        G.add_node(node_index, node=node)  # Attach the node object to the graph

        if parent_index is not None and parent_edge is not None:
            G.add_edge(parent_index, node_index, edge=parent_edge)  # Attach the edge object to the graph
            
            edge_key = parent_edge.reaction_smiles
            if edge_key not in edge_color_mapping:
                edge_color_mapping[edge_key] = len(edge_color_mapping)
            edge_colors[(parent_index, node_index)] = edge_color_mapping[edge_key]

        for edge, child_node in node.subtrees:
            traverse(child_node, node_index, edge, is_from_root=(node_index == 0))

    traverse(root_node)

    fig, ax = plt.subplots(figsize=(15, 10))
    nx.draw(G, pos, labels=labels, with_labels=False, node_size=50, node_color='lightblue', ax=ax, edge_color=[edge_colors[edge] for edge in G.edges()])

    for node_index in pos:
        node = G.nodes[node_index]['node']
        draw_node(ax, node, pos, node_index)

    for edge in G.edges():
        edge_data = G.edges[edge]['edge']
        edge_color = 'C' + str(edge_colors[edge] % 10)  # Cycle through 10 colors
        is_from_root = edge[0] == 0
        draw_edge(ax, edge_data, pos[edge[0]], pos[edge[1]], edge_color, is_from_root)

    plt.savefig('tree.pdf')
    plt.show()

# Example usage
# root = json_to_tree('tree.json')
# visualize_tree(root)
