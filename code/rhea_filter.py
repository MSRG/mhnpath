# Read data from rhea.tsv
with open('rhea.tsv', 'r') as file:
    rhea_data = file.readlines()

# Read data from rhea-obsoletes.tsv
with open('rhea-obsoletes.tsv', 'r') as file:
    obsolete_ids = set(line.strip() for line in file.readlines())

# Filter rows based on obsolete IDs
filtered_rhea_data = [line for line in rhea_data if line.split('\t')[0] not in obsolete_ids]

# Write the filtered data back to rhea.tsv
with open('rhea_filtered.tsv', 'w') as file:
    file.writelines(filtered_rhea_data)

print("Filtered data written to rhea_filtered.tsv")
