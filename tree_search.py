from rdkit import Chem
from rdkit.Chem import AllChem
from price import calculate_cost
from reaction_cond import pred_temperature, pred_solvent_score
from predict import predict
import csv
from mhnreact.inspector import *
import os
import json


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
    dfs_search(start_node, n_enz, n_syn, max_depth, clf_enz, clf_syn1, clf_syn2, clf_syn3, clf_syn4, clf_syn5, start_node, json_pathway)
    print_tree_to_json(start_node, json_pathway)

def dfs_search(node, n_enz, n_syn, max_depth, clf_enz, clf_syn1, clf_syn2, clf_syn3, clf_syn4, clf_syn5, start_node, json_pathway):
    print_tree_to_json(start_node, json_pathway)

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

    node.subtrees.sort(key=lambda x: x[0].score, reverse=True)

    for _, subtree in node.subtrees:
        dfs_search(subtree, max(5, n_enz-10), max(5, n_syn-10), max_depth, clf_enz, clf_syn1, clf_syn2, clf_syn3, clf_syn4, clf_syn5, start_node, json_pathway)

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
                        
find_pathways('Oc(ccc(CC1NCCc(cc2O)c1cc2O)c1)c1O', n_enz=5, n_syn=5, max_depth=5, json_pathway='tree_1.json', device='cpu')

#device can be 'cuda' or 'cpu'
