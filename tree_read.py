import json

class Node:
    def __init__(self, smiles, cost_usd_per_g, depth):
        self.smiles = smiles
        self.cost_usd_per_g = cost_usd_per_g
        self.depth = depth
        self.subtrees = []  # list of tuples (Edge, Node)

class Edge:
    def __init__(self, reaction_smiles, temperature, enzyme, score, rule, label):
        self.reaction_smiles = reaction_smiles
        self.temperature = temperature
        self.enzyme = enzyme
        self.score = score
        self.rule = rule
        self.label = label

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
        root_node = build_node(data)
    except KeyError as e:
        print(f"Error: Missing key '{e}' in the JSON data.")
        return None
    
    return root_node

def build_node(node_data):
    smiles = node_data['smiles']
    cost_usd_per_g = node_data['cost_usd_per_g']
    depth = node_data['depth']
    
    node = Node(smiles, cost_usd_per_g, depth)
    
    if 'subtrees' in node_data and isinstance(node_data['subtrees'], list):
        for subtree_data in node_data['subtrees']:
            edge = build_edge(subtree_data)
            subtree_node = build_node(subtree_data['subtree'])
            node.subtrees.append((edge, subtree_node))
    
    return node

def build_edge(edge_data):
    reaction_smiles = edge_data['reaction_smiles']
    temperature = edge_data['temperature']
    enzyme = edge_data['enzyme']
    score = edge_data['score']
    rule = edge_data['rule']
    label = edge_data['label']
    
    edge = Edge(reaction_smiles, temperature, enzyme, score, rule, label)
    return edge

# Example usage with a JSON file path
# json_file_path = 'tree.json'
# root = json_to_tree(json_file_path)

# visualize_tree(root)
