from mhnreact.inspector import *
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

def predict(smiles, n_enz, n_syn):

    enz_path = 'data/enz_mhn_shuffled.csv'
    syn_path = 'data/syn_small_mhn_shuffled.csv'

    # Load the models
    clf_enz = load_clf('enz_final.pt', model_type='mhn', device='cpu')
    clf_syn = load_clf('syn_finale.pt', model_type='mhn', device='cpu')

    # Predict the reaction
    preds_enz = np.array(clf_enz.forward_smiles([smiles]).tolist()[0])
    preds_syn = np.array(clf_syn.forward_smiles([smiles]).tolist()[0])

    # Get the top n predictions
    top_n_enz = np.argsort(preds_enz)[-n_enz:]
    top_n_syn = np.argsort(preds_syn)[-n_syn:]

    df_enz = pd.read_csv(enz_path, usecols=['label', 'reaction_smarts'])
    df_syn = pd.read_csv(syn_path, usecols=['label', 'reaction_smarts'])

    enz_rules = []
    syn_rules = []

    for label in top_n_enz:
        enz_rules.append(df_enz[df_enz['label'] == label]['reaction_smarts'].values[0])

    for label in top_n_syn:
        syn_rules.append(df_syn[df_syn['label'] == label]['reaction_smarts'].values[0])

    return enz_rules, syn_rules, top_n_enz, top_n_syn

# enz = 'C=C(C)C(CCC(C)=O)CC(=O)O'
# syn = '[CH3:1][O:2][C:3]1[CH:4]=[C:5]2[C:9](=[CH:10][C:11]=1[O:12][CH3:13])[C:8](=[O:14])[CH:21]([C:18]1[CH:17]=[CH:16][N:15]=[CH:20][CH:19]=1)[C:6]2=[CH2:7]'

# r1, r2 = predict(enz)
# print(len(r1))
# print(len(r2))

# r1, r2 = predict(syn)
# print(len(r1))
# print(len(r2))

def evaluate1():

    top1_enz = 0
    top10_enz = 0
    top50_enz = 0
    top100_enz = 0
    top1_syn = 0
    top10_syn = 0
    top50_syn = 0
    top100_syn = 0

    total_enz = 0
    total_syn = 0

    enz_path = 'data/enz_mhn_shuffled.csv'
    syn_path = 'data/syn_small_mhn_shuffled.csv'

    # Load the models
    clf_enz = load_clf('enz_final.pt', model_type='mhn', device='cpu')
    clf_syn = load_clf('syn_finale.pt', model_type='mhn', device='cpu')

    df_enz = pd.read_csv(enz_path, usecols=['label', 'reaction_smarts', 'prod_smiles', 'split'])
    df_syn = pd.read_csv(syn_path, usecols=['label', 'reaction_smarts', 'prod_smiles', 'split'])

    # Filter the DataFrame where split == 'test'
    test_df_enz = df_enz[df_enz['split'] == 'test']
    test_df_syn = df_syn[df_syn['split'] == 'test']

    # Create a dictionary from the filtered DataFrame
    prod_label_dict_enz = dict(zip(test_df_enz['prod_smiles'], test_df_enz['label']))
    prod_label_dict_syn = dict(zip(test_df_syn['prod_smiles'], test_df_syn['label']))

    for prod_smiles, label in prod_label_dict_enz.items():
        preds_enz = np.array(clf_enz.forward_smiles([prod_smiles]).tolist()[0])
        if label in np.flip(np.argsort(preds_enz))[:1]:
            top1_enz += 1
        if label in np.flip(np.argsort(preds_enz))[:10]:
            top10_enz += 1
        if label in np.flip(np.argsort(preds_enz))[:50]:
            top50_enz += 1
        if label in np.flip(np.argsort(preds_enz))[:100]:
            top100_enz += 1
        total_enz += 1

    # print enz metrics
    print('Enzyme metrics:')
    print('Top 1 accuracy: ', top1_enz / total_enz)
    print('Top 10 accuracy: ', top10_enz / total_enz)
    print('Top 50 accuracy: ', top50_enz / total_enz)
    print('Top 100 accuracy: ', top100_enz / total_enz)
    print('starting evaluation, len: ', len(prod_label_dict_syn))
    for prod_smiles, label in prod_label_dict_syn.items():
        preds_syn = np.array(clf_syn.forward_smiles([prod_smiles]).tolist()[0])
        if label in np.flip(np.argsort(preds_syn))[:1]:
            top1_syn += 1
        if label in np.flip(np.argsort(preds_syn))[:10]:
            top10_syn += 1
        if label in np.flip(np.argsort(preds_syn))[:50]:
            top50_syn += 1
        if label in np.flip(np.argsort(preds_syn))[:100]:
            top100_syn += 1
        total_syn += 1
        print(total_syn)

    # print syn metrics
    print('Synthetic metrics:')
    print('Top 1 accuracy: ', top1_syn / total_syn)
    print('Top 10 accuracy: ', top10_syn / total_syn)
    print('Top 50 accuracy: ', top50_syn / total_syn)
    print('Top 100 accuracy: ', top100_syn / total_syn)


def evaluate2():

    top1_enz = []
    top10_enz = []
    top50_enz = []
    top100_enz = []
    top1_syn = []
    top10_syn = []
    top50_syn = []
    top100_syn = []

    total_enz = 0
    total_syn = 0

    enz_path = 'data/enz_mhn_shuffled.csv'
    syn_path = 'data/syn_small_mhn_shuffled.csv'

    # Load the models
    clf_enz = load_clf('enz_final.pt', model_type='mhn', device='cpu')
    clf_syn = load_clf('syn_finale.pt', model_type='mhn', device='cpu')

    df_enz = pd.read_csv(enz_path, usecols=['label', 'reaction_smarts', 'prod_smiles', 'split'])
    df_syn = pd.read_csv(syn_path, usecols=['label', 'reaction_smarts', 'prod_smiles', 'split'])

    # Filter the DataFrame where split == 'test'
    test_df_enz = df_enz[df_enz['split'] == 'test']
    test_df_syn = df_syn[df_syn['split'] == 'test']

    # Create a dictionary from the filtered DataFrame
    prod_label_dict_enz = dict(zip(test_df_enz['prod_smiles'], test_df_enz['label']))
    prod_label_dict_syn = dict(zip(test_df_syn['prod_smiles'], test_df_syn['label']))

    print('starting evaluation, len: ', len(prod_label_dict_enz))

    for prod_smiles, label in prod_label_dict_enz.items():
        preds_enz = np.array(clf_enz.forward_smiles([prod_smiles]).tolist()[0])
        top_enz = np.argsort(preds_enz)[-100:]
        enz_rules = []
        for lab in top_enz:
            enz_rules.append(df_enz[df_enz['label'] == lab]['reaction_smarts'].values[0])
        for rule in enz_rules[:1]:
            out = apply_rule(rule, prod_smiles)
            if out != 0:
                top1_enz.append(1)
            else:
                top1_enz.append(0)
        num10 = 0
        for rule in enz_rules[:10]:
            out = apply_rule(rule, prod_smiles)
            if out != 0:
                num10 += 1
        top10_enz.append(num10)
        num50 = 0
        for rule in enz_rules[:50]:
            out = apply_rule(rule, prod_smiles)
            if out != 0:
                num50 += 1
        top50_enz.append(num50)
        num100 = 0
        for rule in enz_rules[:100]:
            out = apply_rule(rule, prod_smiles)
            if out != 0:
                num100 += 1
        top100_enz.append(num100)
        total_enz += 1
        print(total_enz)

    # print enz metrics
    print('Enzyme metrics:')
    print('Average Number of correct rules:')
    print('Top 1: ', sum(top1_enz) / total_enz)
    print('Top 10: ', sum(top10_enz) / total_enz)
    print('Top 50: ', sum(top50_enz) / total_enz)
    print('Top 100: ', sum(top100_enz) / total_enz)
    print('Accuracy by atleast one:')
    print('Top 1: ', sum(top1_enz) / total_enz)
    print('Top 10: ', sum(1 for num in top10_enz if num > 0) / total_enz)
    print('Top 50: ', sum(1 for num in top50_enz if num > 0) / total_enz)
    print('Top 100: ', sum(1 for num in top100_enz if num > 0) / total_enz)

    for prod_smiles, label in prod_label_dict_syn.items():
        preds_syn = np.array(clf_syn.forward_smiles([prod_smiles]).tolist()[0])
        top_syn = np.argsort(preds_syn)[-100:]
        syn_rules = []
        for lab in top_syn:
            syn_rules.append(df_syn[df_syn['label'] == lab]['reaction_smarts'].values[0])
        for rule in syn_rules[:1]:
            out = apply_rule(rule, prod_smiles)
            if out != 0:
                top1_syn.append(1)
            else:
                top1_syn.append(0)
        num10 = 0
        for rule in syn_rules[:10]:
            out = apply_rule(rule, prod_smiles)
            if out != 0:
                num10 += 1
        top10_syn.append(num10)
        num50 = 0
        for rule in syn_rules[:50]:
            out = apply_rule(rule, prod_smiles)
            if out != 0:
                num50 += 1
        top50_syn.append(num50)
        num100 = 0
        for rule in syn_rules[:100]:
            out = apply_rule(rule, prod_smiles)
            if out != 0:
                num100 += 1
        top100_syn.append(num100)
        total_syn += 1

    # print syn metrics
    print('Synthetic metrics:')
    print('Average Number of correct rules:')
    print('Top 1: ', sum(top1_syn) / total_syn)
    print('Top 10: ', sum(top10_syn) / total_syn)
    print('Top 50: ', sum(top50_syn) / total_syn)
    print('Top 100: ', sum(top100_syn) / total_syn)
    print('Accuracy by atleast one:')
    print('Top 1: ', sum(top1_syn) / total_syn)
    print('Top 10: ', sum(1 for num in top10_syn if num > 0) / total_syn)
    print('Top 50: ', sum(1 for num in top50_syn if num > 0) / total_syn)
    print('Top 100: ', sum(1 for num in top100_syn if num > 0) / total_syn)

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
    
evaluate1()