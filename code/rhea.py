import requests
from bs4 import BeautifulSoup
import pandas as pd

# Load the existing dataset.csv
dataset = pd.read_csv('dataset.csv')

# Load the rhea.tsv file
rhea_data = pd.read_csv('rhea.tsv', sep='\t', header=None, names=['id', 'reaction'])

num = 0

def scrap(id):
    url = "https://www.rhea-db.org/rhea/" + id

    # Send a GET request to the URL
    response = requests.get(url)

    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        # Parse HTML content with BeautifulSoup
        soup = BeautifulSoup(response.text, 'html.parser')

        # Find the heading with text "Enzymes"
        enzymes_heading = soup.find('h4', text='Enzymes')

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
                        
                        # Find the <span> element right next to the <a> element
                        span_element = a_element.find_next_sibling('span')

                        if span_element and "ec" in data_molid:
                            name = span_element.get_text(strip=True)
                            
                            # Return the extracted information
                            return [data_molid[3:], name]
            else:
                return 0
        else:
            return 0
    else:
        print(f"Failed to retrieve data. Status code: {response.status_code}")
        return 0
    return 0

# Iterate over ids in rhea_data
for idx, row in rhea_data.iterrows():
    current_id = row['id']
    num += 1
    # Run scrap(id)
    result = scrap(str(current_id))
    print(num)
    # Process the result
    if result != 0:
        # Append id to id column of dataset.csv
        dataset = dataset._append({'id': current_id, 'source_id': current_id, 'smiles': row['reaction'],
                                  'enzyme_catalysed': 1, 'enzyme_name': result[1], 'ec_number': result[0],
                                  'source': 'RHEA', 'reaction_in_words': ''}, ignore_index=True)
    else:
        # Append id to id column of dataset.csv
        dataset = dataset._append({'id': current_id, 'source_id': current_id, 'smiles': row['reaction'],
                                  'enzyme_catalysed': 0, 'enzyme_name': '', 'ec_number': '',
                                  'source': 'RHEA', 'reaction_in_words': ''}, ignore_index=True)

    # Save the updated dataset.csv
    dataset.to_csv('dataset.csv', index=False)