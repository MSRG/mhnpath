import csv

def compare_and_update_output(output_file, bkms_with_smiles_file):
    # Load the existing data from output.csv
    existing_data = []
    with open(output_file, newline='', encoding='utf-8') as f:
        reader = csv.reader(f)
        next(reader)  # Skip the header
        for row in reader:
            existing_data.append(row)

    # Load data from bkms_with_smiles.txt
    new_data = []
    with open(bkms_with_smiles_file, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            ec_number = row['EC_Number']
            existing_smiles = find_smiles_for_ec(existing_data, ec_number)

            # Compare existing EC numbers with new EC numbers
            if existing_smiles is None:
                new_data.append([len(existing_data) + 1, ec_number, row['smiles']])
                existing_data.append([len(existing_data) + 1, ec_number, row['smiles']])

    # Update the output file
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['ID', 'EC Number', 'SMILES'])
        for item in existing_data:
            writer.writerow(item)

def find_smiles_for_ec(existing_data, ec_number):
    for row in existing_data:
        if row[1] == ec_number:
            return row[2]
    return None

if __name__ == "__main__":
    output_file = 'output.csv'
    bkms_with_smiles_file = 'bkms_w_SMILES.txt'

    compare_and_update_output(output_file, bkms_with_smiles_file)
