import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from price import calculate_cost
from reaction_cond import pred_temperature, pred_solvent_score
from predict import predict
import csv


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
        self.score = 0
        self.rule = rule
        self.label = label

def find_pathways(input_smiles, n_enz=3, n_syn=3, max_depth=3):
    price = get_price(input_smiles)
    if price is None:
        price = 500
    start_node = Node(input_smiles, price, 0)
    dfs_search(start_node, n_enz, n_syn, max_depth)
    print_tree(start_node)

def dfs_search(node, n_enz, n_syn, max_depth):
    if node.cost_usd_per_g <= 100 or node.depth >= max_depth:
        return

    enz_rules, syn_rules, enz_labels, syn_labels = find_applicable_rules(node.smiles, n_enz, n_syn)
    for i in range(len(enz_rules)):
        rule = enz_rules[i]
        label = enz_labels[i]
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

        new_edge = Edge(reaction_smiles, temperature, 0, 1, rule, label)
        for reactant in new_smiles:
            cost_usd_per_g = get_price(reactant)
            if cost_usd_per_g is None:
                cost_usd_per_g = 500
            score = - (temperature/300) - (cost_usd_per_g/500) + solvent_score
            new_edge.score = score
            new_node = Node(reactant, cost_usd_per_g, node.depth + 1)
            node.subtrees.append((new_edge, new_node))

    for i in range(len(syn_rules)):
        rule = syn_rules[i]
        label =syn_labels[i]
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

        new_edge = Edge(reaction_smiles, temperature, 0, 0, rule, label)
        for reactant in new_smiles:
            cost_usd_per_g = get_price(reactant)
            if cost_usd_per_g is None:
                cost_usd_per_g = 500
            score = - (temperature/300) - (cost_usd_per_g/500) + solvent_score
            new_edge.score = score
            new_node = Node(reactant, cost_usd_per_g, node.depth + 1)
            node.subtrees.append((new_edge, new_node))

    node.subtrees.sort(key=lambda x: x[0].score)

    for _, subtree in node.subtrees:
        dfs_search(subtree, n_enz, n_syn, max_depth)

def find_applicable_rules(smiles, n_enz, n_syn):
    enz_rules, syn_rules, enz_labels, syn_labels = predict([smiles], n_enz, n_syn)
    print(smiles+ ' : NUMBER of rules : ' + str(len(enz_rules)) + '----' + str(len(syn_rules)))
    return enz_rules, syn_rules, enz_labels, syn_labels

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
    csv_file = 'buyables.csv'
    with open(csv_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['smiles'] == smiles:
                return float(row['ppg'])
    c = calculate_cost([smiles])
    return min(c)

def print_tree(node, filename='tree.txt', level=0):
    with open(filename, 'a') as file:
        file.write('  ' * level + str(node.__dict__) + '\n')
        for edge, subtree in node.subtrees:
            file.write('  ' * (level + 1) + str(edge.__dict__) + '\n')
            print_tree(subtree, filename, level + 1)

def find_type_disconnect(smiles):
    m = Chem.MolFromSmiles(smiles)
    patt5 = Chem.MolFromSmarts('c1ccccc1OCC')
    patt3 = Chem.MolFromSmarts('c1ccccc1O')
    patt4 = Chem.MolFromSmarts('c1ccccc1OC')
    if m.HasSubstructMatch(patt5):
        return 5
    elif m.HasSubstructMatch(patt3):
        return 3
    elif m.HasSubstructMatch(patt4):
        return 4
    else:
        return 0


# print(apply_rule('[C:5]-[C;H0;D3;+0:4](-[C;D1;H3:6])=[O;H0;D1;+0:7].[C:1]-[C;H0;D3;+0:2](=[O;D1;H0:3])-[O;H1;D1;+0]>>[C:1]-[C;H0;D3;+0:2](=[O;D1;H0:3])-[C;H0;D4;+0:4](-[C:5])(-[C;D1;H3:6])-[OH;D1;+0:7]', 'C=C(C)C(CCC(C)=O)CC(=O)O'))
            
# print(get_price('C=C(C)C1CCC2(C)OC2C1'))
            
# rule = '[C:2]-[C;H0;D3;+0:1](=[O;D1;H0:3])-[NH;D2;+0:4]-[c:5]>>O-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[NH2;D1;+0:4]-[c:5]'
# product = 'CCN(CC)CC(=O)Nc1c(C)cccc1C'
# sm = apply_rule(rule, product)
# print(sm)
# print(get_reactant_smiles(sm))
            
# find_pathways('C1=CC=C2C(=C1)C=CC=C2N')
                        
find_pathways('SC(CCC1=CC=CC=C1)=C', 20, 50, 3)
