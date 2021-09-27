import pandas as pd
import pathlib

def add_features(feature, data, save):
    parcial_feature = data.copy()
    parcial_feature[f"Unique {feature}"] = parcial_feature.index.isin(parcial_feature.drop_duplicates(f"{feature}", keep =False).index)
    parcial_feature.replace({False: 0, True: 1}, inplace=True)
    parcial_feature_2 = parcial_feature.groupby(['ProteinID']).size().sort_values(ascending=False).reset_index(name=f'Total {feature}')
    feature_process = parcial_feature.loc[parcial_feature[f"Unique {feature}"] == 1].drop(columns=[f"Peptide", "Spectrum"])
    feature_process = feature_process.groupby(['ProteinID']).size().reset_index(name=f'Unique {feature}')
    if 'ProteinID' not in save.columns:
        save["ProteinID"] = parcial_feature_2["ProteinID"]
    # Add unique feature column
    save = save.set_index("ProteinID").join(feature_process.set_index("ProteinID")).reset_index()
    # Add total feature column
    save = save.set_index("ProteinID").join(parcial_feature_2.set_index("ProteinID")).reset_index()
    return save


def quantitative_proteomics(path):
    # get all percolator parsed files
    parsed_files = pathlib.Path(path).glob('*_parsed.tsv')

    # Start Empty dataframe to store all _parsed files
    data = pd.DataFrame({'ProteinID':[],\
                    'Peptide':[], \
                    'Spectrum':[]})
    for file in parsed_files:
        df = pd.read_csv(f"{file}", sep='\t', header=0)
        data = data.append(df)

    total = pd.DataFrame({}) 
    total = add_features('Peptide', data, total)
    total = add_features('Spectrum', data, total)
    total = total.fillna(0).astype({"Unique Peptide": int, "Unique Spectrum": int}).sort_values(by='ProteinID', ascending=False)
    total.to_csv("Total_Proteomics_Quantification.tsv", sep="\t", index=False)

quantitative_proteomics("/mnt/c/Users/Pharaohs.son/Desktop/LAB/Ativo/AnnotaPipeline")