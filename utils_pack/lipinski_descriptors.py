from rdkit import Chem
from rdkit.Chem import Lipinski

def compute_lipinski_descriptors(df, smiles_column):
    # Create empty lists to store the computed descriptors
    mw = []
    logp = []
    hbd = []
    hba = []

    # Iterate over the SMILES in the dataframe
    for smiles in df[smiles_column]:
        # Convert the SMILES to a molecule object
        mol = Chem.MolFromSmiles(smiles)

        # Compute the Lipinski descriptors
        mw.append(Lipinski.MolWt(mol))
        logp.append(Lipinski.MolLogP(mol))
        hbd.append(Lipinski.NumHDonors(mol))
        hba.append(Lipinski.NumHAcceptors(mol))

    # Append the computed descriptors to the dataframe
    df['MolecularWeight'] = mw
    df['LogP'] = logp
    df['NumHDonors'] = hbd
    df['NumHAcceptors'] = hba
    
    return df