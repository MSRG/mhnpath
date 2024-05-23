import pandas as pd

# Read the CSV file into a DataFrame
df = pd.read_csv('output.csv')

# Remove the EC Number column
df = df.drop('EC Number', axis=1)

# Keep only unique entries in the SMILES column
df = df.drop_duplicates(subset='SMILES')

# Reset the index and reassign new IDs
df['ID'] = range(1, len(df) + 1)

# Save the modified DataFrame to a new CSV file
df.to_csv('final.csv', index=False)

print("Processing completed. Modified file saved as 'final.csv'")
