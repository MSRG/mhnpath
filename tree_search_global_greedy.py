import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from price import calculate_cost
from reaction_cond import pred_temperature, pred_solvent_score
from predict import predict
import csv
from type_disconnect import find_type_disconnect
import heapq

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

def find_pathways(input_smiles, n_enz=3, n_syn=3, max_depth=3, X=None, Y=None, num=None):
    price = get_price(input_smiles)
    if price is None:
        price = 500
    start_node = Node(input_smiles, price, 0)
    global_greedy_search(start_node, n_enz, n_syn, max_depth, X, Y, num)
    print_tree(start_node)

def global_greedy_search(start_node, n_enz, n_syn, max_depth, X, Y, num):
    pq = []
    heapq.heappush(pq, (-float('inf'), start_node))  # Start node with arbitrarily high score

    while pq:
        _, node = heapq.heappop(pq)
        if node.cost_usd_per_g <= 100 or node.depth >= max_depth:
            continue

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
                type_dis = get_type_disconnect(X, Y, num, node.smiles, reactant)
                score = - min(temperature/300, 1) - min(cost_usd_per_g/500, 1) + solvent_score + type_dis
                new_edge.score = score
                new_node = Node(reactant, cost_usd_per_g, node.depth + 1)
                node.subtrees.append((new_edge, new_node))
                heapq.heappush(pq, (-score, new_node))

        for i in range(len(syn_rules)):
            rule = syn_rules[i]
            label = syn_labels[i]
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
                type_dis = get_type_disconnect(X, Y, num, node.smiles, reactant)
                score = - min(temperature/300, 1) - min(cost_usd_per_g/500, 1) + solvent_score + type_dis
                new_edge.score = score
                new_node = Node(reactant, cost_usd_per_g, node.depth + 1)
                node.subtrees.append((new_edge, new_node))
                heapq.heappush(pq, (-score, new_node))

def find_applicable_rules(smiles, n_enz, n_syn):
    enz_rules, syn_rules, enz_labels, syn_labels = predict([smiles], n_enz, n_syn)
    print(smiles + ' : NUMBER of rules : ' + str(len(enz_rules)) + '----' + str(len(syn_rules)))
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

def get_type_disconnect(X, Y, num, product, reactant):
    type_dis = find_type_disconnect(X, Y, num, product, reactant)
    if type_dis is None:
        return 0
    return type_dis

def print_tree(node, filename='tree_n1_1.txt', level=0):
    with open(filename, 'a') as file:
        file.write('  ' * level + str(node.__dict__) + '\n')
        for edge, subtree in node.subtrees:
            file.write('  ' * (level + 1) + str(edge.__dict__) + '\n')
            print_tree(subtree, filename, level + 1)

# Example usage
find_pathways('CC(C)(C)[C@@H](CS(C)(=O)=O)Nc1nc(-c2c[nH]c3ncc(F)cc23)ncc1F')
