import pandas as pd
from rdkit import Chem

def remove(my_dataframe, smiles_column='smiles') -> pd.DataFrame:
    
    """
    Removes invalid smiles from a given dataframe and returns it.

    Args:
        my_dataframe (pd.DataFrame): DataFrame containing a smiles column with SMILES strings.
        smiles_column (str): Name of the column containing SMILES strings (default: 'smiles').

    Returns:
        pd.DataFrame: DataFrame with invalid smiles removed.
    """
    
    df_smiles = my_dataframe[smiles_column]
    ok_smiles = []
    delete_smiles = []
    i=0

    for ds in df_smiles:
        try:
            cs = Chem.CanonSmiles(ds)
            ok_smiles.append(cs)
            i=i+1
        except:
            print(f"Invalid smiles at line {i} is: {ds}")
            delete_smiles.append(ds)
    
    print(f"There are {len(delete_smiles)} invalid smiles in the dataset")
    
    if len(delete_smiles) > 0:
        print("Returned modified dataframe with invalid smiles removed!")
    
    # Remove the above invalid smiles from the main dataframe
    new_df = my_dataframe[~my_dataframe[smiles_column].isin(delete_smiles)]
    
    return new_df

# Remove rows with smiles containing salts, multiple fragments and transition metals
def contains_transition_metal(smiles) -> bool:
    """Function to check if a SMILES string contains a transition metal


    Args:
        smiles (str): input a SMILES string.

    Returns:
        bool: Returns True if the SMILES string contains a transition metal, False otherwise.
    """
    transition_metals = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                         'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
                         'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg']
    
    mol = Chem.MolFromSmiles(smiles)
    
    for atom in mol.GetAtoms():
        
        if atom.GetSymbol() in transition_metals:
            return True
        
    return False

# Return achiral SMILES
def achirality(smiles) -> str:
    """Function to make SMILES achiral

    Args:
        smiles (str): input a SMILES string

    Returns:
        str: a SMILES string with chiral centres removed
    """
    
    asmiles = Chem.MolFromSmiles(smiles)
    Chem.RemoveStereochemistry(asmiles)
    
    return Chem.MolToSmiles(asmiles)

# Function to filter the dataframe
def filter_smiles(dataframe, smiles_column='smiles') -> pd.DataFrame:
    """
    Perform sanitation and standardization checks on a DataFrame containing SMILES strings.

    Args:
        df (pandas.DataFrame): A DataFrame with a column named 'smiles'.

    Returns:
        pandas.DataFrame: A sanitized and standardized DataFrame.
    """
    
    df = dataframe.copy()
    
    # Remove invalid SMILES
    filtered_df = remove(df, smiles_column='SMILES')
    
    for index, row in filtered_df.iterrows():
        smiles = row[smiles_column]
        
        # Check for multiple fragments
        if '.' in smiles:
            filtered_df = filtered_df.drop(index)
            continue
        
        # Check for transition metals
        if contains_transition_metal(smiles):
            filtered_df = filtered_df.drop(index)
            continue
    
    # removed chirality
    filtered_df['Achiral_SMILES'] = filtered_df['SMILES'].apply(achirality)
    
    # Reset the index after dropping rows
    filtered_df = filtered_df.reset_index(drop=True)
        
    return filtered_df

# Example usage:
# Assuming 'df' is your dataframe with a column named 'smiles'
# df_filtered = filter_smiles(df)
# print(df_filtered)
