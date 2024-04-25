from rdkit import Chem
import pandas as pd

def remove(my_dataframe, smiles_column='smiles'):
    
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