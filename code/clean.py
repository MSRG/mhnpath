import csv

def clean_output_csv(input_file, output_file):
    with open(input_file, newline='') as input_csv, open(output_file, 'w', newline='') as output_csv:
        reader = csv.DictReader(input_csv)
        fieldnames = reader.fieldnames

        # Remove entries with empty "SMILES"
        cleaned_data = [row for row in reader if row['SMILES'].strip()]

        # Renumber the "ID" column
        for i, row in enumerate(cleaned_data, start=1):
            row['ID'] = str(i)

        writer = csv.DictWriter(output_csv, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(cleaned_data)

if __name__ == "__main__":
    input_file = 'output.csv'
    output_file = 'cleaned.csv'

    clean_output_csv(input_file, output_file)
    print("Output CSV cleaned and renumbered.")
