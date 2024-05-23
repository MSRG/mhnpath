import pandas as pd

# Read the original CSV file
original_df = pd.read_csv('/Users/shivesh/Downloads/enzyme/enz_shuffled.csv')

# Randomly select 1000 rows
random_sample_df = original_df.sample(n=300, random_state=42)  # Adjust random_state for reproducibility

# Save the randomly selected rows to a new CSV file
random_sample_df.to_csv('/Users/shivesh/Downloads/enzyme/enz_random.csv', index=False)

# import pandas as pd

# # Read the original CSV file
# original_df = pd.read_csv('/Users/shivesh/Downloads/enzyme/enz_finale.csv')

# # Shuffle the entire DataFrame
# shuffled_df = original_df.sample(frac=1, random_state=42)  # frac=1 shuffles the entire DataFrame

# # Save the shuffled DataFrame to a new CSV file
# shuffled_df.to_csv('/Users/shivesh/Downloads/enzyme/enz_shuffled.csv', index=False)


# import pandas as pd

# # Read the original CSV files
# bkms_df = pd.read_csv('/Users/shivesh/Downloads/enzyme/bkms_finale.csv')
# rhea_df = pd.read_csv('/Users/shivesh/Downloads/enzyme/rhea_finale.csv')

# # Append the data points from rhea to bkms
# enz_df = bkms_df._append(rhea_df, ignore_index=True)

# # Rewrite the 'id' column to start from 1 and increase by 1 for every row
# enz_df['id'] = enz_df.index + 1

# # Save the combined DataFrame to a new CSV file
# enz_df.to_csv('/Users/shivesh/Downloads/enzyme/enz_finale.csv', index=False)

# Path to the CSV file
# path = 'enz_shuffled.csv'

# # Read the CSV file
# df = pd.read_csv(path)

# # Initialize products and rules arrays
# products = []
# rules = []

# # Iterate over each row in the DataFrame
# for index, row in df.iterrows():
#     # Split the 'smiles' column by '>>'
#     smiles_parts = row['smiles'].split('>>')

#     # Extract the part to the right of '>>' and add to the products array
#     product_part = smiles_parts[1].strip()
#     if product_part not in products:
#         products.append(product_part)

#     # Add 'rule_smarts' to the rules array if it doesn't exist
#     if type(row['rule_smarts']) == float:
#         continue
#     rule_smart = row['rule_smarts'].strip()
#     if rule_smart not in rules:
#         rules.append(rule_smart)

# print("Number of unique products:", len(products))
# print("Number of unique rules:", len(rules))

# len_pr = [len(prod) for prod in products]
# len_rl = [len(rule) for rule in rules]

# print("Average length of products:", sum(len_pr) / len(len_pr))
# print("Average length of rules:", sum(len_rl) / len(len_rl))

# products = [prod for prod in products if len(prod) < 600]
# rules = [rule for rule in rules if len(rule) < 800]

# # Print the final products and rules arrays
# print("Number of unique products:", len(products))
# print("Number of unique rules:", len(rules))