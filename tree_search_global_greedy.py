import csv
import json
import heapq
from rdkit import Chem
from predict import predict
from rdkit.Chem import AllChem
from mhnreact.inspector import *
from price import calculate_cost
from reaction_cond import pred_temperature, pred_solvent_score


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

def find_pathways(input_smiles, n_enz=3, n_syn=3, max_depth=3, json_pathway='tree_1.json', device='cpu'):
    # Load the models
    clf_enz = load_clf('enz_final.pt', model_type='mhn', device=device)
    clf_syn1 = load_clf('syn1_final.pt', model_type='mhn', device=device)
    clf_syn2 = load_clf('syn2_final.pt', model_type='mhn', device=device)
    clf_syn3 = load_clf('syn3_final.pt', model_type='mhn', device=device)
    clf_syn4 = load_clf('syn4_final.pt', model_type='mhn', device=device)
    clf_syn5 = load_clf('syn5_final.pt', model_type='mhn', device=device)
    print('getting price')
    price = get_price(input_smiles)
    print('got price: ', price)
    if price is None:
        price = 50000
    start_node = Node(input_smiles, price, 0)
    global_greedy_search(start_node, n_enz, n_syn, max_depth, clf_enz, clf_syn1, clf_syn2, clf_syn3, clf_syn4, clf_syn5, start_node, json_pathway)
    print_tree_to_json(start_node, json_pathway)

def global_greedy_search(node, n_enz, n_syn, max_depth, clf_enz, clf_syn1, clf_syn2, clf_syn3, clf_syn4, clf_syn5, start_node, json_pathway):
    pq = []
    heapq.heappush(pq, (-float('inf'), start_node))  # Start node with arbitrarily high score

    while pq:
        print_tree_to_json(start_node, json_pathway)
        _, node = heapq.heappop(pq)

        if node.cost_usd_per_g <= 100 or node.depth >= max_depth:
            return

        enz_rules, syn_rules, enz_labels = find_applicable_rules(node.smiles, n_enz, n_syn, clf_enz, clf_syn1, clf_syn2, clf_syn3, clf_syn4, clf_syn5)
        for i in range(len(enz_rules)):
            rule = enz_rules[i]
            label = enz_labels[i]
            # label = 0
            reaction_smiles = apply_rule(rule, node.smiles)
            if reaction_smiles == 0:
                continue
            try:
                temperature = pred_temperature(reaction_smiles)
            except:
                temperature = 300
            try:
                solvent_score = pred_solvent_score(reaction_smiles)
            except:
                solvent_score = 0
            new_smiles = get_reactant_smiles(reaction_smiles)

            if start_node.smiles in new_smiles:
                continue

            new_edge = Edge(reaction_smiles, temperature, -1000, 1, rule, label, 0)
            for reactant in new_smiles:
                cost_usd_per_g = get_price(reactant)
                if cost_usd_per_g is None:
                    cost_usd_per_g = 50000
                score = - (temperature/300) - (cost_usd_per_g/500) + solvent_score
                new_edge.score = max(score, new_edge.score)
                new_node = Node(reactant, cost_usd_per_g, node.depth + 1)
                node.subtrees.append((new_edge, new_node))
                heapq.heappush(pq, (-score, new_node))

        for i in range(len(syn_rules)):
            rule = syn_rules[i]
            reaction_smiles = apply_rule(rule, node.smiles)
            if reaction_smiles == 0:
                continue
            try:
                temperature = pred_temperature(reaction_smiles)
            except:
                temperature = 300
            try:
                solvent_score = pred_solvent_score(reaction_smiles)
            except:
                solvent_score = 0
            new_smiles = get_reactant_smiles(reaction_smiles)

            if start_node.smiles in new_smiles:
                continue

            new_edge = Edge(reaction_smiles, temperature, -1000, 0, rule, 0, 0)
            for reactant in new_smiles:
                cost_usd_per_g = get_price(reactant)
                if cost_usd_per_g is None:
                    cost_usd_per_g = 50000
                score = - (temperature/300) - (cost_usd_per_g/500) + solvent_score
                new_edge.score = max(score, new_edge.score)
                new_node = Node(reactant, cost_usd_per_g, node.depth + 1)
                node.subtrees.append((new_edge, new_node))
                heapq.heappush(pq, (-score, new_node))

def find_applicable_rules(smiles, n_enz, n_syn, clf_enz, clf_syn1, clf_syn2, clf_syn3, clf_syn4, clf_syn5):
    enz_rules, syn_rules, enz_labels = predict([smiles], n_enz, n_syn, clf_enz, clf_syn1, clf_syn2, clf_syn3, clf_syn4, clf_syn5)
    print(smiles+ ' : NUMBER of rules : ' + str(len(enz_rules)) + '----' + str(len(syn_rules)))
    return enz_rules, syn_rules, enz_labels

def apply_rule(rule, product):
    rule = '(' + rule.replace('>>', ')>>')
    prod = Chem.MolFromSmiles(product)
    rxn = AllChem.ReactionFromSmarts(rule)
    try:
        reactant = ''
        res = rxn.RunReactants([prod])
        for i in range(len(res[0])):
            if reactant != '':
                reactant += '.'
            reactant += Chem.MolToSmiles(res[0][i])
        return reactant + '>>' + product
    except Exception as e:
        return 0

def get_reactant_smiles(reaction_smiles):
    smiles_parts = reaction_smiles.split('>>')
    reactant_part = smiles_parts[0].strip()
    reactant_parts = reactant_part.split('.')
    return reactant_parts
    
def get_price(smiles):
    try: 
        csv_file = 'buyables.csv'
        compound_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
        with open(csv_file, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                row_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(row['smiles']))
                if row_smiles == compound_smiles:
                    return float(row['ppg'])
    except:
        pass
    csv_file = 'buyables.csv'
    with open(csv_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['smiles'] == smiles:
                return float(row['ppg'])
    c = calculate_cost([smiles])
    return min(c)

def print_tree_to_json(node, filename='tree.json', level=0):
    tree_dict = {
        'smiles': node.smiles,
        'cost_usd_per_g': node.cost_usd_per_g,
        'depth': node.depth,
        'subtrees': []
    }
    
    for edge, subtree in node.subtrees:
        edge_dict = {
            'reaction_smiles': edge.reaction_smiles,
            'temperature': edge.temperature,
            'enzyme': edge.enzyme,
            'score': edge.score,
            'rule': edge.rule,
            'label': edge.label
        }
        subtree_dict = print_tree_to_json(subtree, filename, level + 1)
        edge_dict['subtree'] = subtree_dict
        tree_dict['subtrees'].append(edge_dict)
    
    if level == 0:
        with open(filename, 'w') as file:
            json.dump(tree_dict, file, indent=2)
    else:
        return tree_dict
                        
# Example usage   
# find_pathways('Oc(ccc(CC1NCCc(cc2O)c1cc2O)c1)c1O', n_enz=5, n_syn=5, max_depth=5, json_pathway='tree.json', device='cpu')

# -------------------------
# Command-Line Interface
# -------------------------
#
# This script can also be executed from the command line.
#
# Usage example:
#
#   python tree_search_global_greedy.py -product "Oc(ccc(CC1NCCc(cc2O)c1cc2O)c1)c1O" -n_enz 5 -n_syn 5 -max_depth 5 -json_pathway "tree.json" -device "cpu"
#
# Parameters:
#   -product      : SMILES string of the target product. (Required)
#   -n_enz        : Number of enzyme reaction rules to consider. (Optional, default: 3)
#   -n_syn        : Number of synthetic reaction rules to consider. (Optional, default: 3)
#   -max_depth    : Maximum depth for the tree search. (Optional, default: 3)
#   -json_pathway : Filename to which the resulting pathway tree will be saved in JSON format. (Optional, default: "tree.json")
#   -device       : Device to run the model on; either "cpu" or "cuda". (Optional, default: "cpu")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Tree Search for Pathways: Generate pathways by applying enzyme and synthetic reaction rules."
    )
    parser.add_argument("-product", type=str, required=True,
                        help="SMILES string of the target product (e.g., 'Oc(ccc(CC1NCCc(cc2O)c1cc2O)c1)c1O').")
    parser.add_argument("-n_enz", type=int, default=3,
                        help="Number of enzyme reaction rules to consider (default: 3).")
    parser.add_argument("-n_syn", type=int, default=3,
                        help="Number of synthetic reaction rules to consider (default: 3).")
    parser.add_argument("-max_depth", type=int, default=3,
                        help="Maximum search depth for pathway exploration (default: 3).")
    parser.add_argument("-json_pathway", type=str, default="tree.json",
                        help="Filename for the output JSON file (default: 'tree.json').")
    parser.add_argument("-device", type=str, default="cpu",
                        help="Device to run the model on: 'cpu' or 'cuda' (default: 'cpu').")

    args = parser.parse_args()

    find_pathways(args.product, n_enz=args.n_enz, n_syn=args.n_syn,
                  max_depth=args.max_depth, json_pathway=args.json_pathway, device=args.device)
