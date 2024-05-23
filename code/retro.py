import csv

def append_to_output(output_file, retro_file):
    # Load existing data from output.csv
    existing_data = []
    with open(output_file, newline='', encoding='utf-8') as f:
        reader = csv.reader(f)
        next(reader)  # Skip the header
        for row in reader:
            existing_data.append(row)

    # Load data from retro.csv and append to output.csv
    last_id = int(existing_data[-1][0]) if existing_data else 0
    with open(retro_file, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            last_id += 1
            substrate1_smiles = row['Substrate 1 SMILES']
            substrate2_smiles = row['Substrate 2 SMILES']
            product1_smiles = row['Product 1 SMILES']

            if substrate2_smiles == '':
                smiles = substrate1_smiles + '>>' + product1_smiles
            else:
                smiles = substrate1_smiles + '.' + substrate2_smiles + '>>' + product1_smiles
            
            new_row = [str(last_id), '', smiles]
            existing_data.append(new_row)

    # Update the output file
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['ID', 'EC Number', 'SMILES'])
        for item in existing_data:
            writer.writerow(item)

if __name__ == "__main__":
    output_file = 'output.csv'
    retro_file = 'retro.csv'

    append_to_output(output_file, retro_file)
