#!/usr/bin/python3

# --- A AVENTURA VAI COMEÇAR ---------------------------------------------------

import pandas as pd
from pathlib import Path

# --- SET FILE LOCATION --------------------------------------------------------

base_dir = Path("/mnt/hdd_files/mestrado/pangenoma_leishmania/data")
file_path = Path(base_dir / "13i_OrthoFinder_annotated_proteins/AfterOrthologyProtocol")

# --- OPEN DATAFRAMES ----------------------------------------------------------

orthogroups_df = pd.read_csv(
    Path(file_path / "Orthogroups.tsv"),
    sep="\t"
)

# --- PARSE DATAFRAMES ---------------------------------------------------------

df = pd.read_csv("Total_Proteoma.txt", sep="\t", skiprows=1, index_col=False)
df = df[df["protein"].str.contains("DECOY_")==False]