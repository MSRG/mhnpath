import os
import glob
import pandas as pd
from typing import Optional, List
from rdkit import Chem


# ---------------------------------------------------------------------------
# Path to the folder(s) containing your Mcule CSV files.
# Update MCULE_CSV_DIRS to point at whatever directories you unzip into.
# Every *.csv inside those directories will be searched.
# ---------------------------------------------------------------------------
MCULE_CSV_DIRS = [
    r"C:/Users/dpsso/Downloads/mcule"
]

# Fallback buyables database (used by get_price in tree_search_global_greedy.py)
BUYABLES_CSV = "data/buyables.csv"
 
 
# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------
 
def _canonical(smiles: str) -> Optional[str]:
    """Return RDKit canonical SMILES, or None if the string is unparseable."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol)
    except Exception:
        return None
 
 
def _load_mcule_dataframes() -> pd.DataFrame:
    """
    Read every CSV file found under MCULE_CSV_DIRS and concatenate them.
 
    Expected columns (case-insensitive): 'smiles', 'best_price'
    The 'Mcule ID' column is carried along but not required for lookup.
    """
    frames = []
    for directory in MCULE_CSV_DIRS:
        pattern = os.path.join(directory, "**", "*.csv")
        for csv_path in glob.glob(pattern, recursive=True):
            try:
                df = pd.read_csv(csv_path)
                # Normalise column names to lowercase
                df.columns = [c.strip().lower() for c in df.columns]
                if "smiles" not in df.columns or "best_price" not in df.columns:
                    print(f"[price] Skipping {csv_path}: missing 'smiles' or 'best_price' column")
                    continue
                df = df[["smiles", "best_price"]].dropna(subset=["smiles"])
                frames.append(df)
            except Exception as e:
                print(f"[price] Could not read {csv_path}: {e}")
 
    if not frames:
        return pd.DataFrame(columns=["smiles", "best_price"])
 
    combined = pd.concat(frames, ignore_index=True)
    # Pre-compute canonical SMILES for the whole catalogue once at load time
    combined["canonical"] = combined["smiles"].apply(_canonical)
    return combined
 
 
# Load once at import time so repeated calls are fast
_MCULE_DB: Optional[pd.DataFrame] = None
 
 
def _get_db() -> pd.DataFrame:
    global _MCULE_DB
    if _MCULE_DB is None:
        _MCULE_DB = _load_mcule_dataframes()
        print(f"[price] Loaded {len(_MCULE_DB):,} Mcule entries from CSV files.")
    return _MCULE_DB
 
 
# ---------------------------------------------------------------------------
# Public API  (drop-in replacement for the original price.py)
# ---------------------------------------------------------------------------
 
def calculate_cost(smiles_list: List[str], save_path: str = BUYABLES_CSV) -> list:
    """
    Look up chemical purchase costs from local Mcule CSV files.
 
    Replaces the original API-based implementation.  Prices are looked up by
    canonical SMILES match; if no match is found the compound is returned as
    None (unknown price).
 
    Matched entries are also written back into the buyables.csv database so
    that get_price() in tree_search_global_greedy.py can find them on the
    next call without re-scanning the catalogue.
 
    Parameters
    ----------
    smiles_list : list of str
        SMILES strings to price.
    save_path : str
        Path to the persistent buyables CSV (default: 'data/buyables.csv').
 
    Returns
    -------
    list
        USD/gram prices in the same order as smiles_list; None where unknown.
    """
    db = _get_db()
 
    # Load (or create) the persistent buyables database
    if os.path.exists(save_path):
        buyables = pd.read_csv(save_path)
    else:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        buyables = pd.DataFrame(columns=["smiles", "ppg", "source"])
 
    usd_per_g_prices = []
 
    for smiles in smiles_list:
        canon = _canonical(smiles)
        cost = None
 
        if canon is not None and not db.empty:
            # Match on canonical SMILES (handles tautomers / atom-map differences)
            matches = db[db["canonical"] == canon]
            if matches.empty:
                # Fallback: exact string match on the raw SMILES column
                matches = db[db["smiles"] == smiles]
 
            if not matches.empty:
                # Take the lowest available price across all catalogue files
                best_price = pd.to_numeric(matches["best_price"], errors="coerce").min()
                if pd.notna(best_price):
                    cost = float(best_price)
 
        # ── Update the persistent buyables database ──────────────────────
        if cost is not None:
            if smiles in buyables["smiles"].values:
                existing = buyables.loc[buyables["smiles"] == smiles, "ppg"].values[0]
                if cost < existing:
                    buyables.loc[buyables["smiles"] == smiles, "ppg"] = cost
                    buyables.loc[buyables["smiles"] == smiles, "source"] = "mcule_csv"
            else:
                new_row = pd.DataFrame(
                    {"smiles": [smiles], "ppg": [cost], "source": ["mcule_csv"]}
                )
                buyables = pd.concat([buyables, new_row], ignore_index=True)
 
        usd_per_g_prices.append(cost)
 
    # Persist any new/updated entries
    buyables.to_csv(save_path, index=False)
    return usd_per_g_prices