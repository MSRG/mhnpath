from rdkit import Chem

def find_type_disconnect(X, Y, num, product, reactant):

    if X is None or Y is None:
        return None

    if num == 1:
        core = Y + 'C=C' + X + 'C1=CC=CC=C1'
        type5 = 'C=C' + X + 'C1=CC=CC=C1'
        type4 = 'C' + X + 'C1=CC=CC=C1'
        type3 = X + 'C1=CC=CC=C1'
    elif num == 2:
        core = Y + 'N=C' + X + 'C1=CC=CC=C1'
        type5 = 'N=C' + X + 'C1=CC=CC=C1'
        type4 = 'C' + X + 'C1=CC=CC=C1'
        type3 = X + 'C1=CC=CC=C1'
    elif num == 3:
        core = Y + 'N=N' + X + 'C1=CC=CC=C1'
        type5 = 'N=N' + X + 'C1=CC=CC=C1'
        type4 = 'N' + X + 'C1=CC=CC=C1'
        type3 = X + 'C1=CC=CC=C1'
    elif num == 4:
        core = Y + 'C(=C)C' + X + 'C1=CC=CC=C1'
        type5 = 'C(=C)C' + X + 'C1=CC=CC=C1'
        type4 = 'C' + X + 'C1=CC=CC=C1'
        type3 = X + 'C1=CC=CC=C1'
    elif num == 5:
        core = Y + 'C(=O)C' + X + 'C1=CC=CC=C1'
        type5 = 'C(=O)C' + X + 'C1=CC=CC=C1'
        type4 = 'C' + X + 'C1=CC=CC=C1'
        type3 = X + 'C1=CC=CC=C1'
    elif num == 6:
        core = Y + 'CN=' + X + 'C1=CC=CC=C1'
        type5 = 'CN=' + X + 'C1=CC=CC=C1'
        type4 = 'N=' + X + 'C1=CC=CC=C1'
        type3 = X + 'C1=CC=CC=C1'
    elif num == 7:
        core = Y + 'C2=CC=CN=C2' + X + 'C1=CC=CC=C1'
        type5 = 'C2=CC=CN=C2' + X + 'C1=CC=CC=C1'
        type4 = None
        type3 = X + 'C1=CC=CC=C1'
    elif num == 8:
        core = Y + 'CC=' + X + 'C1=CC=CC=C1'
        type5 = 'CC=' + X + 'C1=CC=CC=C1'
        type4 = 'C=' + X + 'C1=CC=CC=C1'
        type3 = X + 'C1=CC=CC=C1'
    elif num == 9:
        core = Y + 'NN' + X + 'C1=CC=CC=C1'
        type5 = 'NN' + X + 'C1=CC=CC=C1'
        type4 = 'N' + X + 'C1=CC=CC=C1'
        type3 = X + 'C1=CC=CC=C1'
    elif num == 10:
        core = Y + 'SC' + X + 'C1=CC=CC=C1'
        type5 = 'SC' + X + 'C1=CC=CC=C1'
        type4 = 'C' + X + 'C1=CC=CC=C1'
        type3 = X + 'C1=CC=CC=C1'
    elif num == 11:
        core = 'O=P(' + Y + ')(O)C' + X + 'C1=CC=CC=C1'
        type5 = 'O=P(O)C' + X + 'C1=CC=CC=C1'
        type4 = 'C' + X + 'C1=CC=CC=C1'
        type3 = X + 'C1=CC=CC=C1'
    elif num == 12:
        core = Y + 'OC' + X + 'C1=CC=CC=C1'
        type5 = 'OC' + X + 'C1=CC=CC=C1'
        type4 = 'C' + X + 'C1=CC=CC=C1'
        type3 = X + 'C1=CC=CC=C1'
    elif num == 13:
        core = 'O=P(' + Y + ')(O)O' + X + 'C1=CC=CC=C1'
        type5 = 'O=P(O)O' + X + 'C1=CC=CC=C1'
        type4 = 'O' + X + 'C1=CC=CC=C1'
        type3 = X + 'C1=CC=CC=C1'
    elif num == 14:
        core = Y + 'CN' + X + 'C1=CC=CC=C1'
        type5 = 'CN' + X + 'C1=CC=CC=C1'
        type4 = 'N' + X + 'C1=CC=CC=C1'
        type3 = X + 'C1=CC=CC=C1'
    elif num == 15:
        core = Y + 'C(=N)C' + X + 'C1=CC=CC=C1'
        type5 = 'C(=N)C' + X + 'C1=CC=CC=C1'
        type4 = 'C' + X + 'C1=CC=CC=C1'
        type3 = X + 'C1=CC=CC=C1'
    elif num == 16:
        core = Y + 'C2OC2' + X + 'C1=CC=CC=C1'
        type5 = 'C2OC2' + X + 'C1=CC=CC=C1'
        type4 = None
        type3 = X + 'C1=CC=CC=C1'
    elif num == 17:
        core = Y + 'C(=C)C(=O)' + X + 'C1=CC=CC=C1'
        type5 = 'C(=C)C(=O)' + X + 'C1=CC=CC=C1'
        type4 = 'C(=O)' + X + 'C1=CC=CC=C1'
        type3 = X + 'C1=CC=CC=C1'
    elif num == 18:
        core = Y + 'C2=CC=CC=C2' + X + 'C1=CC=CC=C1'
        type5 = 'C2=CC=CC=C2' + X + 'C1=CC=CC=C1'
        type4 = None
        type3 = X + 'C1=CC=CC=C1'
    else:
        core = None
        type5 = None
        type4 = None
        type3 = None

    if core is None:
        return None

    product_molecule = Chem.MolFromSmarts(product)
    reactant_molecule = Chem.MolFromSmarts(reactant)
    core_molecule = Chem.MolFromSmarts(core)

    if product_molecule.HasSubstructMatch(core_molecule) and not reactant_molecule.HasSubstructMatch(core_molecule):
        if type5 is not None:
            if reactant_molecule.HasSubstructMatch(Chem.MolFromSmarts(type5)):
                return 5
        if type4 is not None:
            if reactant_molecule.HasSubstructMatch(Chem.MolFromSmarts(type4)):
                return 4
        if type3 is not None:
            if reactant_molecule.HasSubstructMatch(Chem.MolFromSmarts(type3)):
                return 3
    return None

# X = 'C'
# Y = 'S'
# num = 2

# product = 'SN=CCC1=CC=CC=C1'
# reactant5 = 'N=CCC1=CC=CC=C1'
# reactant4 = 'CCC1=CC=CC=C1'
# reactant3 = 'CC1=CC=CC=C1'
# reactant = 'CCC'

# print(find_type_disconnect(X, Y, num, product, reactant5))