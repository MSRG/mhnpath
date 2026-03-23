import os
import zipfile
import shutil
import urllib.request
import requests
from tqdm import tqdm

# Specify the zip file
article_id = 28673540
zip_file = "28673540.zip"

def get_all_files(article_id):
    all_files = []
    page = 1
    while True:
        url = f"https://api.figshare.com/v2/articles/{article_id}/files?page={page}&page_size=100"
        r = requests.get(url)
        r.raise_for_status()
        batch = r.json()
        if not batch:
            break
        all_files.extend(batch)
        if len(batch) < 100:
            break
        page += 1
    return all_files

print("Fetching file list from figshare API...")
all_files = get_all_files(article_id)
print(f"Found {len(all_files)} files\n")

with zipfile.ZipFile(zip_file, 'w', zipfile.ZIP_STORED) as zf:
    for i, f in enumerate(all_files, 1):
        name = f['name']
        url = f['download_url']
        size = f['size']
        print(f"[{i}/{len(all_files)}] {name} ({size:,} bytes)")
        tmp = name + ".tmp"
        try:
            r = requests.get(url, stream=True, allow_redirects=True)
            r.raise_for_status()
            total = int(r.headers.get('content-length', 0))
            done = 0
            with open(tmp, 'wb') as out:
                for chunk in r.iter_content(chunk_size=65536):
                    out.write(chunk)
                    done += len(chunk)
                    if total:
                        print(f"\r  {done:,}/{total:,} bytes ({done/total*100:.1f}%)", end="", flush=True)
            print()
            zf.write(tmp, name)
            os.remove(tmp)
        except Exception as e:
            if os.path.exists(tmp):
                os.remove(tmp)
            print(f"  ERROR: {e}")

print(f"\nDone! Saved to: {zip_file}")

# Define target directories
csv_dir = "data/"
pickle_h5_dir = "NeuralNet_Cont_Model/"
other_files_dir = "data/model/"

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
