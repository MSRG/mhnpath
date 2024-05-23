import pandas as pd

# Read data from 'bkms_w_SMILES.txt'
data = pd.read_csv('bkms_w_SMILES.txt', sep='\t')

# Create a new DataFrame with the required columns
dataset = pd.DataFrame(columns=["id", "smiles", "enzyme_catalysed", "enzyme_name", "ec_number", "source", "source_id", "reaction_in_words"])

# Copy data to the new DataFrame
dataset["id"] = range(1, len(data) + 1)
dataset["enzyme_catalysed"] = 1
dataset["source"] = "BKMS"
dataset["source_id"] = data["ID"]
dataset["ec_number"] = data["EC_Number"]
dataset["enzyme_name"] = data["Recommended_Name"]
dataset["reaction_in_words"] = data["Reaction"]
dataset["smiles"] = data["smiles"]

# Save the new DataFrame to 'dataset.csv'
dataset.to_csv('dataset.csv', index=False)
