import pandas as pd
from rdkit import Chem

# Remove rows with smiles containing salts, multiple fragments and transition metals

# Function to check if a SMILES string contains a transition metal
def contains_transition_metal(smiles):
    transition_metals = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                         'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
                         'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg']
    
    mol = Chem.MolFromSmiles(smiles)
    
    for atom in mol.GetAtoms():
        
        if atom.GetSymbol() in transition_metals:
            return True
        
    return False

# Function to filter the dataframe
def filter_smiles(dataframe, smiles_column='smiles'):
    """
    Perform sanitation and standardization checks on a DataFrame containing SMILES strings.

    Args:
        df (pandas.DataFrame): A DataFrame with a column named 'smiles'.

    Returns:
        pandas.DataFrame: A sanitized and standardized DataFrame.
    """
    
    filtered_df = dataframe.copy()
    
    for index, row in dataframe.iterrows():
        smiles = row[smiles_column]
        
        # Check for multiple fragments
        if '.' in smiles:
            filtered_df = filtered_df.drop(index)
            continue
        
        # Check for transition metals
        if contains_transition_metal(smiles):
            filtered_df = filtered_df.drop(index)
            continue
        
    return filtered_df

# Example usage:
# Assuming 'df' is your dataframe with a column named 'smiles'
# df_filtered = filter_smiles(df)
# print(df_filtered)