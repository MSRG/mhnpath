import numpy as np
import pandas as pd
from transformers import AutoTokenizer
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern
from sklearn.metrics import mean_squared_error
import joblib
import sys
import warnings

warnings.filterwarnings("error")

# Redirect print statements to a file
sys.stdout = open('enz_gpr_15_2_metrics.txt', 'w')

# Path to the CSV file
path = 'enz_random.csv'

mapped_path = 'enz_random_mapped.csv'

# Read the CSV file
df = pd.read_csv(path)

print("Read CSV")

# Initialize products and rules arrays
products = []
rules = []

# Iterate over each row in the DataFrame
for index, row in df.iterrows():
    # Split the 'smiles' column by '>>'
    smiles_parts = row['smiles'].split('>>')

    if type(row['rule_smarts']) == float:
        continue

    # Extract the part to the right of '>>' and add to the products array
    product_part = smiles_parts[1].strip()
    if product_part not in products:
        products.append(product_part)

    # Add 'rule_smarts' to the rules array if it doesn't exist
    rule_smart = row['rule_smarts'].strip()
    if rule_smart not in rules:
        rules.append(rule_smart)

print("Extracted products and rules")

products = np.array([prod for prod in products if len(prod) < 500])
rules = np.array([rule for rule in rules if len(rule) < 500])

print("length of products: ", len(products))
print("length of rules: ", len(rules))

def encode(string):
    model_name = "DeepChem/ChemBERTa-77M-MLM"
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    return tokenizer.encode(string, padding='max_length', max_length=500)

# Specify the Matern kernel with length_scale and nu parameters
length_scale = 1.0
nu = 1.5
matern_kernel = Matern(length_scale=length_scale, nu=nu)

# Initialize the Gaussian Process model with the Matern kernel
model = GaussianProcessRegressor(kernel=matern_kernel, random_state=42)

batch_size = 20
num_batches = int((0.9 * len(products)) // batch_size)

rule_encoded = [encode(rule) for rule in rules]

print("Encoded rules")

start_idx = len(products) - batch_size
end_idx = len(products)

test_products = products[start_idx:end_idx]
mapped_test_str = pd.read_csv(mapped_path, usecols=[1], skiprows=start_idx, nrows=batch_size).values
mapped_test = np.array([np.array(list(map(int, row[0].split()))) for row in mapped_test_str])
rule_encoded_test = [encode(rule) for rule in rules]
product_encoded_test = [encode(prod) for prod in test_products]

X_test = np.array([[prod_enc] + rule_encoded for prod_enc in product_encoded_test])
Y_test = mapped_test

X_test_flat = X_test.reshape((X_test.shape[0], -1))

print("Processed test data")

for i in range(num_batches):
    if i == 225 or i == 447:
        continue
    start_idx = i * batch_size
    end_idx = (i + 1) * batch_size

    batch_products = products[start_idx:end_idx]

    # Process the batch
    mapped_batch_str = pd.read_csv(mapped_path, usecols=[1], skiprows=start_idx, nrows=batch_size).values
    mapped_batch = np.array([np.array(list(map(int, row[0].split()))) for row in mapped_test_str])
    product_encoded_batch = [encode(prod) for prod in batch_products]

    X_batch = np.array([[prod_enc] + rule_encoded for prod_enc in product_encoded_batch])
    Y_batch = mapped_batch

    # Flatten the encoded input
    x_flat = X_batch.reshape((X_batch.shape[0], -1))

    print(f'Fitting batch {i + 1}')
    try:
        model.fit(x_flat, Y_batch)
        print("fitted: ", i + 1)
    except:
        print('Error fitting batch', i + 1)
        pass

    y_pred, sigma = model.predict(X_test_flat, return_std=True)
    mse = mean_squared_error(Y_test, y_pred)

    print(f'Batch {i + 1} - MSE: {mse}')

    joblib.dump(model, 'enz_gaussian_process_model_15_temp.joblib')
    print("saved temp model")
    del mapped_batch, mapped_batch_str, product_encoded_batch, X_batch, Y_batch, x_flat, batch_products

joblib.dump(model, 'enz_gaussian_process_model_15.joblib')
print("saved model")
print("starting evaluation")
# Define normalization function
def normalize(arr, t_min, t_max):
    for i in range(len(arr)):
        if arr[i] is None:
            arr[i] = 0
    norm_arr = []
    diff = t_max - t_min
    diff_arr = max(arr) - min(arr)
    if diff_arr == 0:
        return arr
    for i in arr:
        temp = (((i - min(arr))*diff)/diff_arr) + t_min
        norm_arr.append(temp)
    return norm_arr

# Define evaluation metrics variables
total_mse = 0
total_t1 = 0
total_t10 = 0
total_num_corr = []

# Loop through batches
for start_idx in range(num_batches * batch_size, len(products), batch_size):
    end_idx = min(start_idx + batch_size, len(products))

    # Get test data for the current batch
    test_products = products[start_idx:end_idx]
    mapped_test_str = pd.read_csv(mapped_path, usecols=[1], skiprows=start_idx, nrows=batch_size).values
    mapped_test = np.array([np.array(list(map(int, row[0].split()))) for row in mapped_test_str])

    product_encoded_test = [encode(prod) for prod in test_products]

    X_test = np.array([[prod_enc] + rule_encoded for prod_enc in product_encoded_test])
    Y_test = mapped_test

    X_test_flat = X_test.reshape((X_test.shape[0], -1))

    # Make predictions on the test set using the model
    y_pred, sigma = model.predict(X_test_flat, return_std=True)

    # Normalize predictions
    y_p = y_pred.copy()
    for i in range(len(y_p)):
        range_to_normalize = (0,1)
        norm = normalize(y_pred[i], range_to_normalize[0], range_to_normalize[1])
        y_p[i] = norm

    # Calculate Mean Squared Error for the batch
    mse = mean_squared_error(Y_test, y_p)
    total_mse += mse

    # Calculate top 1 accuracy, top 10@1 accuracy, and average number of correct rules in top 10 for the batch
    t1 = 0
    t10 = 0
    num_corr = []
    for i in range(len(y_p)):
        sorted_indices = np.argsort(y_p[i])[::-1]

        top_10_indices = sorted_indices[:10]
        top_10_values = y_p[i][top_10_indices]

        for k in range(10):
            if top_10_values[k] < 1:
                break

        if k > 0:
            bl = [Y_test[i][top_10_indices[j]] == 1 for j in range(k)]
            if sum(bl) > 0:
                t1 += 1
        elif Y_test[i][top_10_indices[0]] == 1:
            t1 += 1

        bools = [Y_test[i][j] == 1 for j in top_10_indices]

        num_corr.append(sum(bools))

        if num_corr[-1] > 0:
            t10 += 1

    total_t1 += t1
    total_t10 += t10
    total_num_corr.extend(num_corr)

    del test_products, mapped_test_str, mapped_test, product_encoded_test, X_test, Y_test, X_test_flat, y_pred, sigma, y_p

# Calculate cumulative evaluation metrics over the entire test set
avg_mse = total_mse / (len(products) - num_batches * batch_size)
avg_t1 = total_t1 / (len(products) - num_batches * batch_size)
avg_t10 = total_t10 / (len(products) - num_batches * batch_size)
avg_num_corr = np.mean(total_num_corr)

# Print cumulative evaluation metrics
print('Cumulative Mean Squared Error:', avg_mse)
print('Cumulative Top 1 accuracy:', avg_t1)
print('Cumulative Top 10@1 accuracy:', avg_t10)
print('Cumulative Average number of correct rules in top 10:', avg_num_corr)

# Restore stdout
sys.stdout.close()
sys.stdout = sys.__stdout__
