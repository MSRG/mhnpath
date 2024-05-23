import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

# Path to the CSV file
path = 'enz_shuffled.csv'

out_path = 'enz_shuffled_mapped.csv'

# Read the CSV file
df = pd.read_csv(path)

id = 1

# Iterate over each row in the DataFrame
for index, row in df.iterrows():
    # Split the 'smiles' column by '>>'
    smiles_parts = row['smiles'].split('>>')

    if type(row['rule_smarts']) == float:
        continue

    # Extract the part to the right of '>>' and add to the products array
    product_part = smiles_parts[1].strip()
    if '.' in product_part:
        product_part = product_part.split('.')[0]
    else:
        product = product_part

    reactant = smiles_parts[0].strip()    

    # Add 'rule_smarts' to the rules array if it doesn't exist
    rule = row['rule_smarts'].strip()

    


products = [prod for prod in products if len(prod) < 500]
rules = [rule for rule in rules if len(rule) < 500]

# Open the output CSV file for writing
with open(out_path, 'w') as f:
    f.write("Product,Map\n")  # Write header

    # Iterate over each product and write the row to the CSV file
    for product in products:
        product_map = np.zeros(len(rules), dtype=int)
        try:
            prod = Chem.MolFromSmiles(product)
            j = 0
            for rule in rules:
                rule = '(' + rule.replace('>>', ')>>')
                rxn = AllChem.ReactionFromSmarts(rule)
                try:
                    res = rxn.RunReactants([prod])
                except Exception as e:
                    print(e)
                    res = None
                if res:
                    product_map[j] = 1
                j += 1
        finally:
            # Write the row to the CSV file
            f.write(f"{product},{' '.join(map(str, product_map))}\n")
        