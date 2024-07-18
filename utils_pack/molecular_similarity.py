from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

class similarity:
    
    def __init__(self, mol_1, mol_2):
        self.mol_1 = mol_1
        self.mol_2 = mol_2
        
    def tanimoto_similarity(self) -> float:
        """
        Calculates the Tanimoto similarity between two molecules.
        
        Returns:
            float: Tanimoto similarity between the two molecules.
        """
        
        mol_1 = Chem.MolFromSmiles(self.mol_1)
        mol_2 = Chem.MolFromSmiles(self.mol_2)
        
        fp1 = AllChem.GetMorganFingerprint(mol_1, 2)
        fp2 = AllChem.GetMorganFingerprint(mol_2, 2)
        
        return DataStructs.TanimotoSimilarity(fp1, fp2)
    
    def sample_def(self) -> str:
        return "This is a sample method"