import pandas as pd
import numpy as np
from rdkit import Chem
import splito
from splito import ScaffoldSplit

def train_test_scaffoldsplit(my_dataframe, features, label, smiles_column='SMILES',test_size=None, random_state=None):
    """Generates training and test sets using scaffold split.

    Args:
        my_dataframe (pd.DataFrame): input dataframe.
        features (numpy.ndarray OR pandas.core.series.Series): column/s containing the features.
        label (pandas.core.series.Series): name of the column containing the labels or experimental values.
        smiles_column (str, optional): _description_. Defaults to 'SMILES'.
        test_size (float): proportion of the dataset to include in the test split.
        random_state (int): random seed.

    Returns:
        pd.DataFrame: returns 4 numpy arrays: X_train, X_test, y_train, y_test
    """
    
    # Initialize a splitter
    splitter = ScaffoldSplit(smiles=my_dataframe[smiles_column].tolist(), n_jobs=-1, test_size=test_size, random_state=random_state)

    # Generate indices for training set and test set
    train_idx, test_idx = next(splitter.split(X=my_dataframe[smiles_column].values))

    #assert train_idx.shape[0] > test_idx.shape[0]
    my_dataframe.loc[train_idx, "ScaffoldSplit"] = "train"
    my_dataframe.loc[test_idx, "ScaffoldSplit"] = "test"
    my_dataframe["scaffold"] = splitter.scaffolds
    
    # Split the dataframe based on ScaffoldSplit column
    train_df = my_dataframe[my_dataframe.ScaffoldSplit == "train"]
    test_df = my_dataframe[my_dataframe.ScaffoldSplit == "test"]
    
    X_train = np.array(list(train_df[features]))
    X_test = np.array(list(test_df[features]))
    y_train = train_df[label]
    y_test = test_df[label]
    
    return X_train, X_test, y_train, y_test