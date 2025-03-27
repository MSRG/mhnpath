from chemprice import PriceCollector
import numpy as np
import pandas as pd
import yaml

def calculate_cost(smiles_list, save_path='buyables.csv'):
    # Load API keys from the configuration YAML file
    with open('config.yaml', 'r') as file:
        config = yaml.safe_load(file)

    # Load existing data
    data = pd.read_csv(save_path)
    
    # Initialize the PriceCollector and set API keys
    pc = PriceCollector()
    pc.setMolportApiKey(config.get('molport_api_key'))
    pc.setMCuleApiKey(config.get('mcule_api_key'))
    pc.setChemSpaceApiKey(config.get('chemspace_api_key'))

    # Check status and connection
    print(pc.status())
    print(pc.check())
    
    # Collect prices for the given SMILES list
    all_prices = pc.collect(smiles_list)
    
    # Select the best price
    best = pc.selectBest(all_prices)
    
    # Initialize an array to store USD/g prices
    usd_per_g_prices = []
    
    # Extract USD/g prices
    for smiles in smiles_list:
        price_row = best[best['Input SMILES'] == smiles]
        if len(price_row) == 0 or price_row['USD/g'].values[0] is None:
            cost = None
        else:
            cost = price_row['USD/g'].values[0]
            new_row = pd.DataFrame({'smiles': [smiles], 'ppg': [cost], 'source': ['API']})
            
            # Check if the SMILES already exists in the data
            if smiles in data['smiles'].values:
                # Get the existing row
                existing_row = data[data['smiles'] == smiles]
                existing_price = existing_row['ppg'].values[0]
                
                # Update the row if the new price is lower
                if cost < existing_price:
                    data.loc[data['smiles'] == smiles, 'ppg'] = cost
                    data.loc[data['smiles'] == smiles, 'source'] = 'API'
            else:
                # Append the new row if SMILES does not exist
                data = data.append(new_row, ignore_index=True)
        
        usd_per_g_prices.append(cost)
    
    # Save the updated data back to the CSV file
    data.to_csv(save_path, index=False)
    return usd_per_g_prices
