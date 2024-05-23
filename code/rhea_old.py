import csv

def append_to_output(output_file, rhea_file):
    # Load existing data from output.csv
    existing_data = []
    with open(output_file, newline='', encoding='utf-8') as f:
        reader = csv.reader(f)
        next(reader)  # Skip the header
        for row in reader:
            existing_data.append(row)

    # Load data from rhea.tsv and append to output.csv
    last_id = int(existing_data[-1][0]) if existing_data else 0
    with open(rhea_file, newline='', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            last_id += 1
            new_row = [str(last_id), '', row[1]]
            existing_data.append(new_row)

    # Update the output file
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['ID', 'EC Number', 'SMILES'])
        for item in existing_data:
            writer.writerow(item)

if __name__ == "__main__":
    output_file = 'output.csv'
    rhea_file = 'rhea.tsv'

    append_to_output(output_file, rhea_file)
