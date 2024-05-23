import pandas as pd

# Read the original CSV file
df = pd.read_csv('superNatural.csv', sep=';')

# Create a new 'id' column with numbering starting from 'NP000000001'
df['id'] = ['NP' + str(i).zfill(9) for i in range(1, len(df) + 1)]

# Create a new DataFrame with only 'id' and 'smiles' columns
new_df = df[['id', 'smiles']]

# Save the new DataFrame to a new CSV file
new_df.to_csv('moi.csv', sep=';', index=False)
