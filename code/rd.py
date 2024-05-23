# import rdchiral.main as rd

# # Directly use SMILES/SMARTS
# reaction_smarts = '[C:1][OH:2]>>[C:1][O:2][C]'
# reactant_smiles = 'OCC(=O)OCCCO'
# outcomes = rd.rdchiralRunText(reaction_smarts, reactant_smiles)
# # print(outcomes)

# # Pre-initialize
# rxn = rd.rdchiralReaction(reaction_smarts)
# reactants = rd.rdchiralReactants(reactant_smiles)
# outcomes = rd.rdchiralRun(rxn, reactants)
# # print(outcomes)

# # Get list of atoms that changed as well
# outcomes, mapped_outcomes = rd.rdchiralRun(rxn, reactants, return_mapped=True)
# print(outcomes, mapped_outcomes)

import rdchiral.template_extractor as rd

reaction = {'_id': 0,
  'products': '[CH3:1][c:2]1[cH:3][c:4]([OH:5])[cH:6][c:7]([OH:8])[c:9]1[C:10](=[O:11])[O:17][c:16]1[cH:15][c:14]([CH3:13])[c:21]([C:22](=[O:23])[OH:24])[c:19]([OH:20])[cH:18]1.[OH2:12]',
  'reactants': '[CH3:1][c:2]1[cH:3][c:4]([OH:5])[cH:6][c:7]([OH:8])[c:9]1[C:10](=[O:11])[OH:12].[CH3:13][c:14]1[cH:15][c:16]([OH:17])[cH:18][c:19]([OH:20])[c:21]1[C:22](=[O:23])[OH:24]'}

template = rd.extract_from_reaction(reaction)

print(template['reaction_smarts'])

# import gzip, json, csv

# # with gzip.open('/Users/shivesh/Downloads/enzyme/uspto-reactions-json.tgz') as f:
# #     reactions = json.load(f)

# with gzip.open('/Users/shivesh/Downloads/enzyme/uspto-templates-json.tgz') as f:
#     rules = json.load(f)

# # print(rules[-1])
# # print(rules[0].keys())
    

# # Function to update 'rule_smarts' in the existing 'uspto_final.csv'
# def update_rule_smarts(uspto_data, reaction_id, rule_smarts, i):
#     for j in range(i, len(rules)):
#         row = uspto_data[j]
#         if int(row['id']) == int(reaction_id):
#             row['rule_smarts'] = rule_smarts
#             return

# # Read existing data from uspto_final.csv
# uspto_data = []
# with open('uspto_final.csv', 'r') as csvfile:
#     reader = csv.DictReader(csvfile)
#     for row in reader:
#         uspto_data.append(row)

# # Iterate through rules and update 'rule_smarts'
# for i in range(len(rules)):
#     rule = rules[i]
#     reaction_id = rule.get('reaction_id')
#     if reaction_id:
#         update_rule_smarts(uspto_data, reaction_id, rule.get('reaction_smarts', ''), i)

# # Write updated data to uspto_finale.csv
# with open('uspto_finale.csv', 'w', newline='') as csvfile:
#     fieldnames = ['id', 'smiles', 'enzyme_catalysed', 'enzyme_name', 'ec_number', 'source', 'source_id', 'reaction_in_words', 'rule_smarts']
#     writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

#     writer.writeheader()
#     writer.writerows(uspto_data)


# Function to create 'uspto_final.csv' based on 'reactions'
# def create_uspto_final_csv(reactions):
#     uspto_data = []

#     for reaction in reactions:
#         uspto_row = {
#             'id': reaction['_id'],
#             'smiles': f"{reaction['reactants']}>>{reaction['products']}",
#             'enzyme_catalysed': 0,
#             'enzyme_name': '',
#             'ec_number': '',
#             'source': 'uspto',
#             'source_id': reaction['source_id'],
#             'reaction_in_words': '',
#             'rule_smarts': '',
#         }

#         uspto_data.append(uspto_row)

#     # Write data to uspto_final.csv
#     with open('uspto_final.csv', 'w', newline='') as csvfile:
#         fieldnames = ['id', 'smiles', 'enzyme_catalysed', 'enzyme_name', 'ec_number', 'source', 'source_id', 'reaction_in_words', 'rule_smarts']
#         writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

#         writer.writeheader()
#         writer.writerows(uspto_data)

# # Call the function with the provided 'reactions' data
# create_uspto_final_csv(reactions)

# # Read existing data from bkms_finale.csv
# existing_data = []
# with open('bkms_finale.csv', 'r') as csvfile:
#     reader = csv.DictReader(csvfile)
#     for row in reader:
#         existing_data.append(row)

# # Function to update 'rule_smarts' in the existing data
# def update_rule_smarts(existing_data, rule_id, rule_smarts):
#     for row in existing_data:
#         if int(row['id']) == int(rule_id):
#             row['rule_smarts'] = rule_smarts
#             print(row)

# # Iterate through rules and update 'rule_smarts'
# for rule in rules:
#     if rule.get('template_set') == 'bkms':
#         for reference_id in rule.get('references', []):
#             update_rule_smarts(existing_data, reference_id, rule['reaction_smarts'])

# # Write updated data back to bkms_finale.csv
# with open('bkms_final.csv', 'w', newline='') as csvfile:
#     fieldnames = ['id', 'smiles', 'enzyme_catalysed', 'enzyme_name', 'ec_number', 'source', 'source_id', 'reaction_in_words', 'rule_smarts']
#     writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

#     writer.writeheader()
#     writer.writerows(existing_data)

# Function to write data to CSV file
# def write_to_csv(filename, data):
#     with open(filename, 'w', newline='') as csvfile:
#         fieldnames = ['id', 'smiles', 'enzyme_catalysed', 'enzyme_name', 'ec_number', 'source', 'source_id', 'reaction_in_words', 'rule_smarts']
#         writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

#         writer.writeheader()

#         for reaction in data:
#             if reaction.get('template_set') == 'bkms':
#                 writer.writerow({
#                     'id': reaction['reaction_id'],
#                     'smiles': reaction['smiles'],
#                     'enzyme_catalysed': '1',
#                     'enzyme_name': reaction['Recommended_Name'],
#                     'ec_number': reaction['EC_Number'],
#                     'source': 'bkms',
#                     'source_id': reaction['reaction_id'],
#                     'reaction_in_words': reaction['Reaction'],
#                     'rule_smarts': '',
#                 })

# # Write data to CSV file
# write_to_csv('bkms_finale.csv', reactions)
    
# def write_to_csv(filename, data):
#     with open(filename, 'w', newline='') as csvfile:
#         fieldnames = ['id', 'smiles', 'enzyme_catalysed', 'enzyme_name', 'ec_number', 'source', 'source_id', 'reaction_in_words', 'rule_smarts']
#         writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

#         writer.writeheader()

#         for rule in data:
#             if rule.get('template_set') == 'reaxys':
#                 writer.writerow({
#                     'id': rule['_id'],
#                     'smiles': '',
#                     'enzyme_catalysed': '0',
#                     'enzyme_name': '',
#                     'ec_number': '',
#                     'source': 'reaxys',
#                     'source_id': rule['_id'],
#                     'reaction_in_words': '',
#                     'rule_smarts': rule['reaction_smarts'],
#                 })

# # Write data to CSV file
# write_to_csv('reaxys_finale.csv', rules)

from rxnmapper import RXNMapper
rxnmapper = RXNMapper()

ex = "Cc1cc(OC(=O)c2c(C)cc(O)cc2O)cc(O)c1C(=O)O.O>>Cc1cc(O)cc(O)c1C(=O)O.Cc1cc(O)cc(O)c1C(=O)O"

res = rxnmapper.get_attention_guided_atom_maps([ex])
res = res[0]['mapped_rxn']
res = res.split('>>')

print(res)
