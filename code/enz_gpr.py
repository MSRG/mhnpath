from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd
from transformers import AutoTokenizer
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, Matern
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
import joblib
import sys
from rdkit import RDLogger

# Redirect print statements to a file
sys.stdout = open('metrics.txt', 'w')

### Section 1: Read the CSV file and extract unique products and rules ###

# Path to the CSV file
path = 'enz_random.csv'

# Read the CSV file
df = pd.read_csv(path)

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

products = [prod for prod in products if len(prod) < 500]
rules = [rule for rule in rules if len(rule) < 500]

# Print the final products and rules arrays
print("Number of unique products:", len(products))
print("Number of unique rules:", len(rules))

### Section 2: Create a map of products and rules and make input-output arrays ###

RDLogger.DisableLog('rdApp.*')

print('Making map of products and rules...')
def make_map(products, rules):
    map = np.zeros((len(products),len(rules)), dtype=int)
    i = 0

    for product in products:
        prod = Chem.MolFromSmiles(product)
        j = 0
        for rule in rules:
            rule = '(' + rule.replace('>>', ')>>')
            rxn = AllChem.ReactionFromSmarts(rule)
            try:
                res = rxn.RunReactants([prod])
            except Exception as e:
                print(e)
                res = None
            if res:
                map[i][j] = 1
            j+=1
        i+=1
        print('Mapping: ', product)
    return map

def encode(string):
    model_name = "DeepChem/ChemBERTa-77M-MLM"
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    return tokenizer.encode(string, padding='max_length', max_length=500)

mapped = make_map(products, rules)

rule_encoded = [encode(rule) for rule in rules]

product_encoded = [[encode(prod)] for prod in products]
X = []
for enc in product_encoded:
    X.append(enc + rule_encoded)

Y = mapped

x = np.array(X)
y = np.array(Y)

print("Shape of input array: ", x.shape)
print("Shape of output array: ", y.shape)

### Section 3: Train a Gaussian Process model on the input-output arrays ##

# Flatten the encoded input
x_flat = x.reshape((x.shape[0], -1))

# Split the data into training and testing sets
x_train, x_test, y_train, y_test = train_test_split(x_flat, y, test_size=0.2, random_state=42)

print('Training Gaussian Process model...')
# # Define the Gaussian Process model with an RBF kernel
# kernel = 1.0 * RBF(length_scale=1.0)
# model = GaussianProcessRegressor(kernel=kernel, random_state=42)

# # Train the Gaussian Process model
# model.fit(x_train, y_train)

# joblib.dump(model, 'gaussian_process_model.joblib')

# print('Gaussian Process model trained successfully!')
# ### Section 4: Evaluate the model ###

# # Evaluate the model

# print('Evaluating metrics for Gaussian Process model')
# # Make predictions on the test set
# y_pred, sigma = model.predict(x_test, return_std=True)
# 'k1__constant_value': [1e-08, 1e-06, 1e-05, 1e-04]
# 1.0 * RBF(length_scale=1.0), 
# Define a wider range for the lower bound of k1__constant_value
param_grid = {
    'kernel': [1.0 * RBF(length_scale=1.0), 1.0 * Matern(length_scale=1.0)],
    'alpha': [1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4],
    'k1__constant_value': [1e-08, 1e-06, 1e-05, 1e-04]}

model = GaussianProcessRegressor(random_state=42)

# Define the number of folds for cross-validation
k_folds = 5
kf = KFold(n_splits=k_folds, shuffle=True, random_state=42)

# Perform GridSearchCV with k-fold cross-validation for hyperparameter tuning
grid_search = GridSearchCV(model, param_grid, scoring='neg_mean_squared_error', cv=kf)
grid_search.fit(x_train, y_train)

# Get the best model from the grid search
best_model = grid_search.best_estimator_

# Save the best model to a file
joblib.dump(best_model, 'best_gaussian_process_model.joblib')

# Print the best hyperparameters
print('Best hyperparameters:', grid_search.best_params_)

# Section 4: Evaluate the model

print('Evaluating metrics for the best Gaussian Process model')

# Make predictions on the test set using the best model
y_pred, sigma = best_model.predict(x_test, return_std=True)

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

y_p = y_pred.copy()

for i in range(len(y_p)):
    range_to_normalize = (0,1)
    norm = normalize(y_pred[i], range_to_normalize[0], range_to_normalize[1])
    y_p[i] = norm

print(f'Mean Squared Error: {mean_squared_error(y_test, y_p)}')
        
t10 = 0
num_corr = []
t1 = 0
for i in range(len(y_p)):
    sorted_indices = np.argsort(y_p[i])[::-1]

    top_10_indices = sorted_indices[:10]
    top_10_values = y_p[i][top_10_indices]

    for k in range(10):
        if top_10_values[k] < 1:
            break

    if k > 0:
        bl = [y_test[i][top_10_indices[j]] == 1 for j in range(k)]
        if sum(bl) > 0:
            t1 += 1
    elif y_test[i][top_10_indices[0]] == 1:
        t1 += 1
    
    bools = [y_test[i][j] == 1 for j in top_10_indices]

    num_corr.append(sum(bools))

    if num_corr[-1] > 0:
        t10 += 1

print('Top 1 accuracy:', t1/len(y_p))
print('Top 10@1 accuracy:', t10/len(y_p))
print('Average number of correct rules in top 10:', np.mean(num_corr))

# Restore stdout
sys.stdout.close()
sys.stdout = sys.__stdout__
