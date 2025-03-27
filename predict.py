from mhnreact.inspector import *
import numpy as np
import pandas as pd

def predict(smiles, n_enz, n_syn, clf_enz, clf_syn1, clf_syn2, clf_syn3, clf_syn4, clf_syn5):

    enz_path = 'data/enz_mhn_shuffled.csv'
    syn1_path = 'data/syn_mhn_split_1.csv'
    syn2_path = 'data/syn_mhn_split_2.csv'
    syn3_path = 'data/syn_mhn_split_3.csv'
    syn4_path = 'data/syn_mhn_split_4.csv'
    syn5_path = 'data/syn_mhn_split_5.csv'


    # # Load the models
    # clf_enz = load_clf('enz_final.pt', model_type='mhn', device='cpu')
    # clf_syn1 = load_clf('syn1_final.pt', model_type='mhn', device='cpu')
    # clf_syn2 = load_clf('syn2_final.pt', model_type='mhn', device='cpu')
    # clf_syn3 = load_clf('syn3_final.pt', model_type='mhn', device='cpu')
    # clf_syn4 = load_clf('syn4_final.pt', model_type='mhn', device='cpu')
    # clf_syn5 = load_clf('syn5_final.pt', model_type='mhn', device='cpu')

    # Predict the reaction
    preds_enz = np.array(clf_enz.forward_smiles([smiles]).tolist()[0])
    preds_syn1 = np.array(clf_syn1.forward_smiles([smiles]).tolist()[0])
    preds_syn2 = np.array(clf_syn2.forward_smiles([smiles]).tolist()[0])
    preds_syn3 = np.array(clf_syn3.forward_smiles([smiles]).tolist()[0])
    preds_syn4 = np.array(clf_syn4.forward_smiles([smiles]).tolist()[0])
    preds_syn5 = np.array(clf_syn5.forward_smiles([smiles]).tolist()[0])

    # Get the top n predictions
    top_n_enz = np.flip(np.argsort(preds_enz))[:n_enz]
    top_n_syn1 = np.flip(np.argsort(preds_syn1))[:n_syn]
    top_n_syn2 = np.flip(np.argsort(preds_syn2))[:n_syn]
    top_n_syn3 = np.flip(np.argsort(preds_syn3))[:n_syn]
    top_n_syn4 = np.flip(np.argsort(preds_syn4))[:n_syn]
    top_n_syn5 = np.flip(np.argsort(preds_syn5))[:n_syn]

    df_enz = pd.read_csv(enz_path, usecols=['label', 'reaction_smarts'])
    df_syn1 = pd.read_csv(syn1_path, usecols=['label', 'reaction_smarts'])
    df_syn2 = pd.read_csv(syn2_path, usecols=['label', 'reaction_smarts'])
    df_syn3 = pd.read_csv(syn3_path, usecols=['label', 'reaction_smarts'])
    df_syn4 = pd.read_csv(syn4_path, usecols=['label', 'reaction_smarts'])
    df_syn5 = pd.read_csv(syn5_path, usecols=['label', 'reaction_smarts'])

    enz_rules = []
    syn_rules = {}

    for label in top_n_enz:
        enz_rules.append(df_enz[df_enz['label'] == label]['reaction_smarts'].values[0])

    for i in range(len(top_n_syn1)):
        label = top_n_syn1[i]
        rule = df_syn1[df_syn1['label'] == label]['reaction_smarts'].values[0]
        score = preds_syn1[label]
        if rule in syn_rules:
            syn_rules[rule] = max(syn_rules[rule], score)
        else:
            syn_rules[rule] = score

        label = top_n_syn2[i]
        rule = df_syn2[df_syn2['label'] == label]['reaction_smarts'].values[0]
        score = preds_syn2[label]
        if rule in syn_rules:
            syn_rules[rule] = max(syn_rules[rule], score)
        else:
            syn_rules[rule] = score

        label = top_n_syn3[i]
        rule = df_syn3[df_syn3['label'] == label]['reaction_smarts'].values[0]
        score = preds_syn3[label]
        if rule in syn_rules:
            syn_rules[rule] = max(syn_rules[rule], score)
        else:
            syn_rules[rule] = score

        label = top_n_syn4[i]
        rule = df_syn4[df_syn4['label'] == label]['reaction_smarts'].values[0]
        score = preds_syn4[label]
        if rule in syn_rules:
            syn_rules[rule] = max(syn_rules[rule], score)
        else:
            syn_rules[rule] = score

        label = top_n_syn5[i]
        rule = df_syn5[df_syn5['label'] == label]['reaction_smarts'].values[0]
        score = preds_syn5[label]
        if rule in syn_rules:
            syn_rules[rule] = max(syn_rules[rule], score)
        else:
            syn_rules[rule] = score

    sorted_syn_rules = sorted(syn_rules, key=syn_rules.get, reverse=True)

    return enz_rules, sorted_syn_rules[:n_syn], top_n_enz
