import gzip, json, csv
with gzip.open('/Users/shivesh/Downloads/enzyme/bkms-and-reaxys-templates.json.gz') as f:
    reactions = json.load(f)

print(list(reactions[171000].keys()))
print(reactions[171002])

# Define the output CSV file name
# output_csv_file = "rules_uspto.csv"

# # Write data to CSV file
# with open(output_csv_file, mode='w', newline='') as csv_file:
#     fieldnames = ['id', 'rule_smarts', 'source']
#     writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

#     # Write header
#     writer.writeheader()

#     # Write data
#     for idx, entry in enumerate(reactions, start=1):
#         if 'reaction_smarts' in entry:
#             writer.writerow({'id': idx, 'rule_smarts': entry['reaction_smarts'], 'source': 'uspto'})