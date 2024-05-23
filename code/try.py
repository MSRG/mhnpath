# import pandas as pd

# path = 'syn_shuffled.csv'

# df = pd.read_csv(path)

# # Group by the "rule_smarts" column and count occurrences
# rule_counts = df.groupby('rule_smarts').size()

# # Filter for rules that occur more than 10 times
# frequent_rules = rule_counts[rule_counts > 5]

# # Get the rules as a DataFrame
# frequent_rules_df = pd.DataFrame({'rule_smarts': frequent_rules.index})

# # Save the frequent rules to a new CSV file
# frequent_rules_df.to_csv('syn_common_rules.csv', index=False)


# import gzip, json, csv

# out = 'syn_freq_rules.csv'

# with gzip.open('/Users/shivesh/Downloads/enzyme/bkms-and-reaxys-templates.json.gz') as f:
#     rules = json.load(f)

# # Write rules with count >= 10 to the CSV file
# with open(out, 'w', newline='') as csvfile:
#     writer = csv.writer(csvfile)
#     writer.writerow(['rules_smarts'])
#     for r in rules:
#         if r['count'] >= 10:
#             writer.writerow([r['reaction_smarts']])

import pandas as pd
import csv

input_file = 'uspto_finale.csv'
output_file = 'syn_freq_rules.csv'

# Read the input CSV file
df = pd.read_csv(input_file)

# Group by the "rule_smarts" column and count occurrences
rule_counts = df.groupby('rule_smarts').size()

# Filter for rules that occur more than 10 times
frequent_rules = rule_counts[rule_counts > 2]

# Check if the output file already contains the frequent rules
existing_rules = set()
try:
    with open(output_file, 'r') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip the header row
        for row in reader:
            existing_rules.add(row[0])
except FileNotFoundError:
    pass

# Append new frequent rules to the output file
with open(output_file, 'a', newline='') as csvfile:
    writer = csv.writer(csvfile)
    for rule in frequent_rules.index:
        if rule not in existing_rules:
            writer.writerow([rule])
