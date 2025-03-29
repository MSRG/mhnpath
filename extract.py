import os
import zipfile
import shutil

# Specify the zip file
zip_file = "28673540.zip"

# Define target directories
csv_dir = "data/"
pickle_h5_dir = "NeuralNet_Cont_Model/"
other_files_dir = "data/models/"

# Ensure all target directories exist
os.makedirs(csv_dir, exist_ok=True)
os.makedirs(pickle_h5_dir, exist_ok=True)
os.makedirs(other_files_dir, exist_ok=True)

# Extract the zip file
with zipfile.ZipFile(zip_file, 'r') as zip_ref:
    extracted_dir = "extracted_files"
    os.makedirs(extracted_dir, exist_ok=True)
    zip_ref.extractall(extracted_dir)
    print(f"Extracted files to {extracted_dir}")

# Move files to their respective directories based on their extensions
for root, _, files in os.walk(extracted_dir):
    for file in files:
        source_path = os.path.join(root, file)

        # Determine the target directory based on file extension
        if file.endswith('.csv'):
            destination_dir = csv_dir
        elif file.endswith('.pickle') or file.endswith('.h5'):
            destination_dir = pickle_h5_dir
        else:
            destination_dir = other_files_dir

        # Ensure the destination directory exists (optional, already created above)
        os.makedirs(destination_dir, exist_ok=True)

        # Move the file
        destination_path = os.path.join(destination_dir, file)
        shutil.move(source_path, destination_path)
        print(f"Moved {file} to {destination_path}")

# Clean up: remove the extracted folder after moving files
shutil.rmtree(extracted_dir)
print(f"Cleaned up extracted directory: {extracted_dir}")

# Delete the zip file after processing
if os.path.exists(zip_file):
    os.remove(zip_file)
    print(f"Deleted zip file: {zip_file}")