# import pandas as pd

# # Read the existing CSV file (moi.csv) with semicolon delimiter
# existing_df = pd.read_csv('moi.csv', sep=';')

# # Read the data from the original TXT file (chembl.txt) with tab separation
# new_data = pd.read_csv('chembl.txt', sep='\t')

# # Create a new 'id' column with numbering starting from 'MD000000001'
# new_data['id'] = ['MD' + str(i).zfill(9) for i in range(1, len(new_data) + 1)]

# # Select only the 'id' and 'canonical_smiles' columns from the new data
# new_data = new_data[['id', 'canonical_smiles']]

# # Rename the 'canonical_smiles' column to 'smiles' in the new data
# new_data = new_data.rename(columns={'canonical_smiles': 'smiles'})

# # Concatenate the existing DataFrame and the new data
# combined_df = pd.concat([existing_df, new_data])

# # Save the combined DataFrame to the same CSV file (moi.csv) with semicolon delimiter
# combined_df.to_csv('moi2.csv', sep=';', index=False)

# from picamera2 import Picamera2
# from libcamera import controls

# picam2.camera_controls['LensPosition']

# import json
# import pandas as pd

# def read_json_to_dict(json_file):
#     with open(json_file, 'r') as f:
#         data_dict = json.load(f)
#     return data_dict

# def append_to_solvent(toxin, solvents_file):
#     data = pd.read_csv(solvents_file)
#     for i in range(len(toxin)):
#         new_row = pd.DataFrame({'smiles': [toxin[i]['moldb_smiles']], 'score': -1})
#         data = data.append(new_row, ignore_index=True)
#     data.to_csv(solvents_file, index=False)

# # Example usage:
# toxins_file = 'toxins.json'
# toxins_dict = read_json_to_dict(toxins_file)
# solvents_file = 'solvent.csv'
# append_to_solvent(toxins_dict, solvents_file)

# import json

# def read_json_to_dict(json_file):
#     with open(json_file, 'r') as f:
#         data_dict = json.load(f)
#     return data_dict

# # Example usage:
# toxins_file = 'toxins.json'
# toxins_dict = read_json_to_dict(toxins_file)
# print(toxins_dict[0]['moldb_smiles'])

# import csv

# def append_to_solvent(csv_file, output_file):
#     with open(csv_file, 'r') as f:
#         reader = csv.DictReader(f, delimiter=';')
#         with open(output_file, 'a', newline='') as outfile:
#             writer = csv.writer(outfile)
#             for row in reader:
#                 if row['smiles'] and row['smiles'] != 'NULL':
#                     writer.writerow([row['smiles'], '1'])

# # Example usage:
# superNatural_file = 'superNatural.csv'
# solvent_file = 'solvent.csv'
# append_to_solvent(superNatural_file, solvent_file)

import csv
import re
from rdkit import Chem
from rdkit.Chem import inchi

def inchikey_to_smiles(inchikey):
    # Convert InChIKey to InChI using PubChem
    from urllib.request import urlopen
    import json

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/InChI/JSON"
    response = urlopen(url)
    data = json.load(response)

    if 'PropertyTable' not in data or 'Properties' not in data['PropertyTable'] or len(data['PropertyTable']['Properties']) == 0:
        raise ValueError(f"Could not find InChI for InChIKey {inchikey}")

    inchi_string = data['PropertyTable']['Properties'][0]['InChI']

    # Convert InChI to SMILES using RDKit
    mol = inchi.MolFromInchi(inchi_string)
    if mol is None:
        raise ValueError(f"Could not convert InChI {inchi_string} to molecule")

    smiles = Chem.MolToSmiles(mol)
    return smiles

# Example usage
# inchikey = 'FTOVXSOBNPWTSH-UHFFFAOYSA-N'
# smiles = inchikey_to_smiles(inchikey)
# print(f'SMILES: {smiles}')

def process_toxins(input_csv, output_csv):
    # Open the input and output CSV files
    with open(input_csv, mode='r', newline='', encoding='utf-8') as infile, \
         open(output_csv, mode='a', newline='', encoding='utf-8') as outfile:
        
        reader = csv.reader(infile)
        writer = csv.writer(outfile)
        
        # Skip the header of the input file
        next(reader)
        
        for row in reader:
            row_text = ','.join(row)
            match = re.search(r'InChIKey=([A-Z\-]+)', row_text)
            
            if match:
                inchikey = match.group(1)
                try:
                    smiles = inchikey_to_smiles(inchikey)
                    writer.writerow([smiles, -1])
                except Exception as e:
                    print(f"Error processing InChIKey {inchikey}: {e}")

# Example usage
input_csv = 'toxins.csv'
output_csv = 'toxicity.csv'
process_toxins(input_csv, output_csv)