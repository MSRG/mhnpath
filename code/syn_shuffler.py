# import pandas as pd

# # Path to the input CSV file
# input_file = 'syn_finale.csv'

# # Path to the output CSV file
# output_file = 'syn_finale.csv'

# # Read the CSV file
# df = pd.read_csv(input_file)

# # Reset the 'id' column to sequential integers starting from 1
# df['id'] = range(1, len(df) + 1)

# # Write the updated DataFrame to a new CSV file
# df.to_csv(output_file, index=False)

# print(f'Updated data saved to {output_file}')

import pandas as pd

# Path to the input CSV file
input_file = 'syn_finale.csv'

# Path to the output CSV file
output_file = 'syn_shuffled.csv'

# Read the CSV file
df = pd.read_csv(input_file)

# Shuffle the rows
shuffled_df = df.sample(frac=1, random_state=42)  # Set a seed (random_state) for reproducibility

# Reset the 'id' column to sequential integers starting from 1
shuffled_df['id'] = range(1, len(shuffled_df) + 1)

# Write the shuffled DataFrame to a new CSV file
shuffled_df.to_csv(output_file, index=False)

print(f'Shuffled data saved to {output_file}')

# import pandas as pd

# # Load data from CSV files
# uspto_data = pd.read_csv('uspto_finale.csv')
# reaxys_data = pd.read_csv('reaxys_finale.csv')

# # Append data from reaxys to uspto
# syn_data = pd.concat([uspto_data, reaxys_data], ignore_index=True)

# # Save the combined data to syn_finale.csv
# syn_data.to_csv('syn_finale.csv', index=False)
