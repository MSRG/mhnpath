import csv
import re

def extract_smiles(smiles_file):
    smiles_dict = {}
    with open(smiles_file, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                compound, smiles = parts
                smiles_dict[compound] = smiles
    return smiles_dict

def process_reactions(input_file, output_file, smiles_dict):
    compound_count = {}
    error_reactions = []
    no_smiles_data = []
    output_data = []  # Define the output_data list

    with open(input_file, newline='', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)  # Skip the header

        current_id = 1
        for row in reader:
            try:
                ec_number = row[1]
                reaction = row[3]

                if ' <=> ' in reaction:
                    reactants, products = reaction.split(' <=> ')
                    p2, r2 = reaction.split(' <=> ')

                    # Split reactants and products into compounds
                    rc = r2.split(' + ')
                    pc = p2.split(' + ')
                    rc3 = []
                    pc3 = []

                    for c in rc:
                        cin = c.split(' ')
                        for item in cin:
                            if item == 'n' or re.match(r'^\d+$', item):
                                cin.remove(item)
                        rc3 += [' '.join(cin)]
                        
                    for c in pc:
                        cin = c.split(' ')
                        for item in cin:
                            if item == 'n' or re.match(r'^\d+$', item):
                                cin.remove(item)
                        pc3 += [' '.join(cin)]

                else:
                    reactants, products = reaction.split(' = ')

                # Split reactants and products into compounds
                rc = reactants.split(' + ')
                pc = products.split(' + ')
                rc2 = []
                pc2 = []

                for c in rc:
                    cin = c.split(' ')
                    for item in cin:
                        if item == 'n' or re.match(r'^\d+$', item):
                            cin.remove(item)
                    rc2 += [' '.join(cin)]
                    
                for c in pc:
                    cin = c.split(' ')
                    for item in cin:
                        if item == 'n' or re.match(r'^\d+$', item):
                            cin.remove(item)
                    pc2 += [' '.join(cin)]

                no_smile = 0

                for compound in rc2 + pc2:
                    if compound not in smiles_dict:
                        no_smile = 1
                        no_smiles_data.append([compound])
                        if compound not in compound_count:
                            compound_count[compound] = 1
                        else:
                            compound_count[compound] += 1

                if no_smile == 0:
                    # Generate the output data
                    smiles = '.'.join([smiles_dict.get(compound, compound) for compound in rc2]) + '>>' + '.'.join([smiles_dict.get(compound, compound) for compound in pc2])
                    output_data.append([current_id, ec_number, smiles])
                    current_id += 1
                    if ' <=> ' in reaction:
                        smiles = '.'.join([smiles_dict.get(compound, compound) for compound in rc3]) + '>>' + '.'.join([smiles_dict.get(compound, compound) for compound in pc3])
                        output_data.append([current_id, ec_number, smiles])
                        current_id += 1
            except Exception as e:
                error_reactions.append([reaction])

    # Write the error reactions to a separate CSV file
    with open('error_reactions.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Reaction'])
        writer.writerows(error_reactions)

    # Write the compounds without SMILES to a separate CSV file
    with open('no_smiles.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Compound', 'Count'])
        for compound, count in compound_count.items():
            if compound not in smiles_dict:
                writer.writerow([compound, count])

    # Save the main data to the output_file
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['ID', 'EC Number', 'SMILES'])
        writer.writerows(output_data)

if __name__ == "__main__":
    ext_names_to_smiles_file = 'ext_names_to_smiles.txt'
    reactions_bkms_file = 'Reactions_BKMS.csv'
    output_file = 'output.csv'

    smiles_dict = extract_smiles(ext_names_to_smiles_file)
    process_reactions(reactions_bkms_file, output_file, smiles_dict)
