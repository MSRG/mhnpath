The data collection phase involves aggregating diverse datasets from BKMS, RHEA, RetroBioCat, USPTO, and other relevant sources to create a comprehensive repository of reaction smiles and chemical information. This data will include enzymatic reactions, standard chemical reactions, and crucial molecules with associated pricing details.

#### To-Do List:
1. Retrieve datasets from BKMS, RHEA, RetroBioCat, and USPTO databases.
   - **Target**: Collect around 50,000 enzymatics reaction smiles and a separate larget dataset from USPTO.
2. Compile reaction smiles and associated data, ensuring integrity and consistency.
   - **Target**: Clean and process data to remove duplicates or inconsistencies.
3. Curate an adversarial dataset comprising enzymatic reactions that do not occur.
   - **Target**: Create a dataset of around 3 per enzymatic reaction, the approach is similar to [this paper](https://www.nature.com/articles/s41467-023-38347-2#Sec14).
4. Create a lookup dataset housing important molecules and their pricing information.
   - **Target**: Collect data on molecules with associated pricing, use ASKOS and MCurate.
