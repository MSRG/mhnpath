import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
import sys

RDLogger.DisableLog('rdApp.*')

# Redirect print statements to a file
sys.stdout = open('syn_mapper_metrics.txt', 'w')
sys.stderr = open('syn_mapper_errors.txt', 'w')

# Path to the CSV file
path = 'syn_shuffled.csv'

out_path = 'syn_shuffled_mapped.csv'

# Read the CSV file
df = pd.read_csv(path)

print('read csv')

# Initialize products and rules arrays
products = []
rules = []

# Iterate over each row in the DataFrame
for index, row in df.iterrows():
    if type(row['rule_smarts']) != str:
        continue

    if type(row['smiles']) == str and row['smiles'] != '':

        # Split the 'smiles' column by '>>'
        smiles_parts = row['smiles'].split('>>')
        # Extract the part to the right of '>>' and add to the products array
        product_part = smiles_parts[1].strip()
        if product_part not in products and len(product_part) < 500:
            products.append(product_part)

    # Add 'rule_smarts' to the rules array if it doesn't exist
    rule_smart = row['rule_smarts'].strip()
    if rule_smart not in rules and len(rule_smart) < 500:
        rules.append(rule_smart)

print(len(products), len(rules))

# Open the output CSV file for writing
with open(out_path, 'w') as f:
    f.write("Product,Map\n")  # Write header

    # Iterate over each product and write the row to the CSV file
    for product in products:
        product_map = np.zeros(len(rules), dtype=int)

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

        # Write the row to the CSV file
        f.write(f"{product},{' '.join(map(str, product_map))}\n")
        print('mapped: ', product)