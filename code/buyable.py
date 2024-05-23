import json
import csv

# Define the JSON file path
json_file_path = 'buyables.json'

# Read JSON data from file
with open(json_file_path, 'r') as jsonfile:
    json_data = jsonfile.read()

data = json.loads(json_data)

# Define the CSV file path
csv_file_path = 'buyables.csv'

# Write the JSON data to CSV
with open(csv_file_path, 'w', newline='') as csvfile:
    fieldnames = ['smiles', 'ppg', 'source']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for item in data:
        writer.writerow(item)

print("CSV file has been created successfully.")
