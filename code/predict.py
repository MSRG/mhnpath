from mhnreact.inspector import *
import numpy as np
import pandas as pd

def predict(smiles, n_enz=3, n_syn=3):

    enz_path = 'data/enz_mhn_shuffled.csv'
    syn_path = 'data/syn_small_mhn_shuffled.csv'

    # Load the models
    clf_enz = load_clf('enz_final.pt', model_type='mhn', device='cpu')
    clf_syn = load_clf('syn_finale.pt', model_type='mhn', device='cpu')

    # Predict the reaction
    preds_enz = np.array(clf_enz.forward_smiles([smiles]).tolist()[0])
    preds_syn = np.array(clf_syn.forward_smiles([smiles]).tolist()[0])

    print(preds_enz)
    print(np.argsort(preds_enz))

    # Get the top n predictions
    top_n_enz = np.flip(np.argsort(preds_enz))[:n_enz]
    top_n_syn = np.flip(np.argsort(preds_syn))[:n_syn]

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

# r1, r2, g1, g2 = predict(enz)
# print(r1)
# print(r2)

# r1, r2 = predict(syn)
# print(len(r1))
# print(len(r2))