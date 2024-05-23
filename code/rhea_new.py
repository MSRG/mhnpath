import requests
from bs4 import BeautifulSoup
import pandas as pd
from rxnmapper import RXNMapper
import rdchiral.template_extractor as rd


rxnmapper = RXNMapper()

# Load the existing dataset.csv
dataset = pd.read_csv('rhea_new.csv')

# Load the rhea.tsv file
rhea_data = pd.read_csv('rhea_filtered.tsv', sep='\t', header=None, names=['id', 'reaction'])

num = 0

def go_ec(id):
    url = "https://www.rhea-db.org/rhea/" + id

    go = 0
    ec = 0

    # Send a GET request to the URL
    response = requests.get(url)

    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        soup = BeautifulSoup(response.text, 'html.parser')
        # Print the text content of the webpage
        # print(response.text)
        # Find the heading with text "Enzymes"
        enzymes_heading = soup.find('h2', string='Cross-references')

        if enzymes_heading:
            # Find the next table after the "Enzymes" heading
            enzymes_table = enzymes_heading.find_next('table')

            if enzymes_table:
                # Find all rows in the table
                rows = enzymes_table.find_all('tr')

                # Extract data-molid and name from each row
                for row in rows:
                    # Find the <a> element within the row
                    a_element = row.find('a', class_='molName')

                    if a_element:
                        data_molid = a_element.get('data-molid')

                        if "GO" in data_molid:
                            go = data_molid
                        
                        if "ec" in data_molid:
                            ec = data_molid[3:]
                return [go, ec]  
        return 0      
    else:
        print(f"Failed to retrieve data. Status code: {response.status_code}")
        return 100

def from_go(numbers, id):
    if numbers == 0:
        return 0
    if numbers == 100:
        return 100

    go = numbers[0]
    ec = numbers[1]

    try:
        if go:
            url = 'https://amigo.geneontology.org/amigo/term/' + go
            response = requests.get(url)

            if response.status_code == 200:
                soup = BeautifulSoup(response.text, 'html.parser')
                name = soup.find('title')
                first = str(name).index('"')
                last = str(name).index(' activity')
                return str(name)[first+1:last]
    except:
        return 0
    return 0

def from_ec(numbers, id):
    if numbers == 0:
        return 0
    if numbers == 100:
        return 100

    go = numbers[0]
    ec = numbers[1]

    if ec:
        url = 'https://enzyme.expasy.org/EC/' + ec
        response = requests.get(url)

        if response.status_code == 200:
            soup = BeautifulSoup(response.text, 'html.parser')

            name = soup.find('strong')
            return name.contents[0][:-1]
    return 0

def scrap(id):
    numbers = go_ec(id)
    if numbers == 0:
        return 0
    if numbers == 100:
        return 'pass'
    name = from_go(numbers, id)
    if name == 0:
        name = from_ec(numbers, id)
    if name == 0:
        return 0
    return [numbers[1], name]

# Iterate over ids in rhea_data
for idx, row in rhea_data.iterrows():
    current_id = row['id']
    num += 1

    # Run scrap(id)
    result = scrap(str(current_id))
    if result == 'pass':
        pass
    print(num)
    smiles = row['reaction']
    try:
        if len(smiles) < 512:
            mapped = rxnmapper.get_attention_guided_atom_maps([smiles])
            mapped = mapped[0]['mapped_rxn']
            mapped = mapped.split('>>')
            reaction = {'_id': 0, 'products': mapped[1], 'reactants': mapped[0]}
            template = rd.extract_from_reaction(reaction)
            template = template['reaction_smarts']
        else:
            template = ''
    except:
        template = ''
    # Process the result
    if result != 0:
        # Append id to id column of dataset.csv
        dataset = dataset._append({'id': current_id, 'source_id': current_id, 'smiles': smiles,
                                  'enzyme_catalysed': 1, 'enzyme_name': result[1], 'ec_number': result[0],
                                  'source': 'RHEA', 'reaction_in_words': '', 'rule_smarts': template}, ignore_index=True)
    else:
        # Append id to id column of dataset.csv
        dataset = dataset._append({'id': current_id, 'source_id': current_id, 'smiles': smiles,
                                  'enzyme_catalysed': 0, 'enzyme_name': '', 'ec_number': '',
                                  'source': 'RHEA', 'reaction_in_words': '', 'rule_smarts': template}, ignore_index=True)

    # Save the updated dataset.csv
    dataset.to_csv('rhea_new.csv', index=False)
