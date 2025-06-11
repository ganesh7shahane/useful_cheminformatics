import streamlit as st

import streamlit as st
import pandas as pd

st.set_page_config(page_title="Molecule Data Analyzer")

uploaded_file = st.file_uploader("Upload your data file", 
                                type=["csv", "xlsx"],
                                help="Supports CSV and Excel files")

if uploaded_file is not None:
    # Read file based on extension
    if uploaded_file.name.endswith('.xlsx'):
        df = pd.read_excel(uploaded_file)
    else:
        df = pd.read_csv(uploaded_file)
    
    # Display basic stats
    st.subheader("Dataset Statistics")
    col1, col2 = st.columns(2)
    col1.metric("Total Rows", df.shape[0])
    col2.metric("Total Columns", df.shape[1])
    
    # Display SMILES molecules
    st.subheader("Top Molecules")
    if 'smiles' in df.columns:
        st.dataframe(df['smiles'].head(10), 
                    use_container_width=True,
                    hide_index=True,
                    column_config={"smiles": "SMILES Structure"})
    else:
        st.error("No 'smiles' column found in the uploaded file")
