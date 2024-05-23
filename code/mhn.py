from mhnreact.inspector import *

list_models()

clf = load_clf(list_models()[12], model_type='mhn', device='cpu')

from mhnreact.data import load_USPTO
X,y = load_USPTO('sm')

preds = clf.forward_smiles(X['test'])

print(X['test'])