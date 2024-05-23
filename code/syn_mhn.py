import pandas as pd
import numpy as np

path = 'syn_shuffled.csv'
out_path = 'syn_small_mhn_shuffled.csv'
freq_path = 'syn_freq_rules.csv'

df = pd.read_csv(path)
df2 = pd.read_csv(freq_path)

r = []
for index, row in df2.iterrows():
    r.append(row['rules_smarts'])

print('length: ', len(r))

# Initialize lists to store data for the new DataFrame
ids = []
prod_smiles = []
reactants_can = []
reaction_smarts = []
split = []
enzyme_catalysed = []
enzyme_name = []
ec_number = []
source = []
source_id = []
reaction_in_words = []
label = []

lc = 0
rc = {}
# Iterate over each row in the DataFrame
for index, row in df.iterrows():
    if len(r) == 0:
        continue
    if type(row['rule_smarts']) == float or type(row['smiles']) == float:
        continue
    if row['rule_smarts'] not in r:
        continue
    r.remove(row['rule_smarts'])
    # Split the 'smiles' column by '>>'
    smiles_parts = row['smiles'].split('>>')

    # Extract the product part
    product_part = smiles_parts[1].strip()
    if '.' in product_part:
        product_part = product_part.split('.')[0]

    reactant = smiles_parts[0].strip()    

    # Add data to lists
    ids.append(index + 1)  # Assuming index starts from 0
    prod_smiles.append(product_part)
    reactants_can.append(reactant)

    rule = row['rule_smarts'].strip()
    if rule in rc:
        label.append(rc[rule])
    else:
        rc[rule] = lc
        label.append(lc)
        lc += 1
    reaction_smarts.append(rule)
    
    # Determine the split value
    split_value = np.random.choice(['train', 'valid', 'test'], p=[0.8, 0.1, 0.1])
    split.append(split_value)
    
    # Add data from original DataFrame
    enzyme_catalysed.append(row['enzyme_catalysed'])
    enzyme_name.append(row['enzyme_name'])
    ec_number.append(row['ec_number'])
    source.append(row['source'])
    source_id.append(row['source_id'])
    reaction_in_words.append(row['reaction_in_words'])

    print('processed: ', index)


# Create a new DataFrame
new_df = pd.DataFrame({
    'id': ids,
    'prod_smiles': prod_smiles,
    'reactants_can': reactants_can,
    'reaction_smarts': reaction_smarts,
    'label': label,
    'split': split,
    'enzyme_catalysed': enzyme_catalysed,
    'enzyme_name': enzyme_name,
    'ec_number': ec_number,
    'source': source,
    'source_id': source_id,
    'reaction_in_words': reaction_in_words
})

# Write the new DataFrame to CSV
new_df.to_csv(out_path, index=False)
