from chemprice import PriceCollector
import numpy as np
import pandas as pd
# import subprocess
# import pricePrediction

def calculate_cost(smiles_list, save_path='buyables.csv'):
    # Load existing data
    data = pd.read_csv(save_path)
    
    # Initialize the PriceCollector and set API keys
    pc = PriceCollector()
    pc.setMolportApiKey("880d8343-8ui2-418c-9g7a-68b4e2e78c8b")
    pc.setMCuleApiKey("75f3cd5b8fdd67dd0f62557bfbddcb62976ab98b")
    pc.setChemSpaceApiKey("88Js8CeWIumkRqRA08xtLBzLYbc_cmy0S8d_8-_q65-JbP6-E1Or_AwIJFf0Lix6")

    # Check status and connection
    print(pc.status())
    print(pc.check())
    
    # Collect prices for the given SMILES list
    all_prices = pc.collect(smiles_list)
    # print(all_prices.to_string())
    
    # Select the best price
    best = pc.selectBest(all_prices)
    # print(best.to_string())
    
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


# def calculate_cost(smiles_list, save_path='buyables.csv'):
#     data = pd.read_csv(save_path)
#     pc = PriceCollector()
#     pc.setMolportApiKey("880d8343-8ui2-418c-9g7a-68b4e2e78c8b")
#     pc.setMCuleApiKey("75f3cd5b8fdd67dd0f62557bfbddcb62976ab98b")
#     # pc.setChemSpaceApiKey("88Js8CeWIumkRqRA08xtLBzLYbc_cmy0S8d_8-_q65-JbP6-E1Or_AwIJFf0Lix6")

#     # Check status and connection
#     print(pc.status())
#     print(pc.check())
    
#     # Collect prices for the given SMILES list
#     all_prices = pc.collect(smiles_list)
#     print(all_prices.to_string())
    
#     # Select the best price
#     best = pc.selectBest(all_prices)
#     print(best.to_string())
    
#     # Initialize an array to store USD/g prices
#     usd_per_g_prices = []
    
#     # Extract USD/g prices
#     for smiles in smiles_list:
#         price_row = best[best['Input SMILES'] == smiles]
#         if len(price_row) == 0 or price_row['USD/g'].values[0] is None:
#             cost = None
#         else:
#             cost = price_row['USD/g'].values[0]
#             new_row = pd.DataFrame({'smiles': [smiles], 'ppg': [cost], 'source': ['API']})
#             data = data.append(new_row, ignore_index=True)
#         usd_per_g_prices.append(cost)
#     data.to_csv(save_path, index=False)
#     return usd_per_g_prices


def predict_cost(path):
    # Construct the command to run the script
    command = ["python", "-m", "pricePrediction.predict", path, "-o", path[:-4] + "_cost.csv"]

    # Call the script using subprocess
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print("Error:", e)

    df = pd.read_csv(path[:-4] + "_cost.csv")
    smiles_list = df['CoPriNet']
    return smiles_list

# path = '/Users/shivesh/Downloads/enzyme/data/testData/set.csv'
# # print(predict_cost(path))
# smiles = ['C1=CC=C2C(=C1)C=CC=C2N', 'Nc1cccc2ccc(C(=O)O)cc12']
# smiles = ['24890knlejqf', '23423']
# smiles = ['O=C(NS(C1=CC=CC=C1[N+]([O-])=O)(=O)=O)C2=CC=CC=C2']
# smiles = ['O[Mn](=O)(=O)=O']
# smiles = ['CCOC(=O)C(C)Br', 'Cc1ccccc1', 'CCO']
# print(calculate_cost(smiles))