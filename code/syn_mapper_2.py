import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
import sys

# Check if the correct number of command-line arguments is provided
if len(sys.argv) != 2:
    print("Usage: python program.py <input_argument>")
    sys.exit(1)

# Extract the input argument from the command-line arguments
path = sys.argv[1]

RDLogger.DisableLog('rdApp.*')

# Redirect print statements to a file
sys.stdout = open(path[:-4] + '_metrics.txt', 'w')
sys.stderr = open(path[:-4] +  '_errors.txt', 'w')

freq_path = 'syn_freq_rules.csv'

out_path = path[:-4] + '_mapped.csv'

# Read the CSV file
df = pd.read_csv(path)

df2 = pd.read_csv(freq_path)

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

    # # Add 'rule_smarts' to the rules array if it doesn't exist
    # rule_smart = row['rule_smarts'].strip()
    # if rule_smart not in rules and len(rule_smart) < 500:
    #     rules.append(rule_smart)
            
for index, row in df2.iterrows():
    rule_smart = row['rules_smarts'].strip()
    if rule_smart not in rules and len(rule_smart) < 500:
        rules.append(rule_smart)

print(len(products), len(rules))

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
        except Exception as e:
            pass

        # Write the row to the CSV file
        f.write(f"{product},{' '.join(map(str, product_map))}\n")
        print('mapped: ', product)
