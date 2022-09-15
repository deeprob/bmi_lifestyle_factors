import os
import numpy as np
import pandas as pd
import json
import subprocess
from sklearn.decomposition import PCA

CURRENT_DIR_PATH = os.path.dirname(os.path.abspath("__file__"))

##########################
# field encodings parser #
##########################

def get_pheno_encoding_filepath(root_dir, pheno_type, pheno_cat):
    """
    This function accepts
    1) root_dir: where the sample to phenotype value table is stored under type -> category -> id hierarchy
    2) pheno_type: the type which the field belongs to
    3) pheno_cat: the category which the field belongs to 
    It returns
    the filepath of the json file that contains fields encodings of all fields 
    that fall under the type and category specified
    """
    pheno_json_path = os.path.join(
        root_dir, pheno_type, pheno_cat, f"fields_data_coding.json"
        )
    assert os.path.exists(pheno_json_path)
    return pheno_json_path


def read_pheno_encodings(pheno_json_path, pheno_field_id):
    pheno_field_id = str(pheno_field_id)
    with open(pheno_json_path, "r") as f:
        field_encoding_dict = json.load(f)
    return field_encoding_dict[pheno_field_id]

##############################
# read the lifestyle factors #
##############################

def get_continuous_factors(lifestyle_df, store_dir, fe_dir):
    lifestyle_df = lifestyle_df.loc[(lifestyle_df.shortlist=="X") & (lifestyle_df.Type=="continuous")]
    dfs = [pd.read_csv(os.path.join(store_dir, t, g, "tables", f"{i}_quantile.csv"), usecols=["eid", "merged"], index_col=["eid"]).rename(columns={"merged": i}) for t,g,i in zip(lifestyle_df.Type, lifestyle_df.Phenotype_group, lifestyle_df.Phenotype_ID)]
    fe_dicts_dict = {i:read_pheno_encodings(get_pheno_encoding_filepath(fe_dir, t, g), i) for t,g,i in zip(lifestyle_df.Type, lifestyle_df.Phenotype_group, lifestyle_df.Phenotype_ID)}
    df = pd.concat(dfs, axis=1).sort_values("eid")
    return df, fe_dicts_dict


def get_integer_factors(lifestyle_df, store_dir, fe_dir):
    lifestyle_df = lifestyle_df.loc[(lifestyle_df.shortlist=="X") & (lifestyle_df.Type=="integer")]
    dfs = [pd.read_csv(os.path.join(store_dir, t, g, "tables", f"{i}_quantile.csv"), usecols=["eid", "merged"], index_col=["eid"]).rename(columns={"merged": i}) for t,g,i in zip(lifestyle_df.Type, lifestyle_df.Phenotype_group, lifestyle_df.Phenotype_ID)]
    fe_dicts_dict = {i:read_pheno_encodings(get_pheno_encoding_filepath(fe_dir, t, g), i) for t,g,i in zip(lifestyle_df.Type, lifestyle_df.Phenotype_group, lifestyle_df.Phenotype_ID)}
    df = pd.concat(dfs, axis=1).sort_values("eid")
    return df, fe_dicts_dict

def get_catsingle_factors(lifestyle_df, store_dir, fe_dir):
    lifestyle_df = lifestyle_df.loc[(lifestyle_df.shortlist=="X") & (lifestyle_df.Type=="categorical_single")]
    dfs = [pd.read_csv(os.path.join(store_dir, t, g, "tables", f"{i}.csv"), usecols=["eid", "merged"], index_col=["eid"]).rename(columns={"merged": i}) for t,g,i in zip(lifestyle_df.Type, lifestyle_df.Phenotype_group, lifestyle_df.Phenotype_ID)]
    fe_dicts_dict = {i:read_pheno_encodings(get_pheno_encoding_filepath(fe_dir, t, g), i) for t,g,i in zip(lifestyle_df.Type, lifestyle_df.Phenotype_group, lifestyle_df.Phenotype_ID)}
    df = pd.concat(dfs, axis=1).sort_values("eid")
    return df, fe_dicts_dict

def get_catmultiple_factors(lifestyle_df, store_dir, fe_dir):
    lifestyle_df = lifestyle_df.loc[(lifestyle_df.shortlist=="X") & (lifestyle_df.Type=="categorical_multiple")]
    dfs = [pd.read_csv(os.path.join(store_dir, t, g, "tables", f"{i}.csv"), usecols=lambda x: x.startswith("eid") or x.startswith("binarized"), index_col=["eid"]).rename(columns=lambda col: col.replace("binarized", f"{i}")) for t,g,i in zip(lifestyle_df.Type, lifestyle_df.Phenotype_group, lifestyle_df.Phenotype_ID)]
    fe_dicts_dict = {i:read_pheno_encodings(get_pheno_encoding_filepath(fe_dir, t, g), i) for t,g,i in zip(lifestyle_df.Type, lifestyle_df.Phenotype_group, lifestyle_df.Phenotype_ID)}
    df = pd.concat(dfs, axis=1).sort_values("eid")
    return df, fe_dicts_dict


#############
# normalize #
#############

def standardize_df(df):
    # df = (df - df.mean())/df.std()
    df = (df-df.min())/(df.max()-df.min())
    return df

def replace_negs(df):
    return df.mask(df<0)

def feature_pca(df):
    pca = PCA(n_components=2)
    Xt = pca.fit_transform(df.values)
    plot_df = pd.DataFrame()
    plot_df.index = df.index
    plot_df["component_1"] = Xt[:,0]
    plot_df["component_2"] = Xt[:,1]
    plot_df = plot_df.merge(df, left_index=True, right_index=True)
    return plot_df

def normalize_single_field_function(standardized_df):
    pca_df = pd.DataFrame()
    pca_df["component_1"] = standardized_df.iloc[:, 0]
    pca_df = pca_df.merge(standardized_df, left_index=True, right_index=True)
    return pca_df

def normalize_integers(int_df, sel_fields):
    # select relevant fields
    selected_df = int_df.loc[:, sel_fields]
    # get rid of negative values
    selected_noneg_df = replace_negs(selected_df).dropna(axis=0)
    # standardize df
    selected_noneg_standardized_df = standardize_df(selected_noneg_df)
    # get pca components
    if len(sel_fields) == 1:
        pca_df = normalize_single_field_function(selected_noneg_standardized_df)
    else:
        pca_df = feature_pca(selected_noneg_standardized_df)    
    return pca_df

def normalize_catsingle_helper(df, field_encodings, onesided_fields):
    # TODO: high and low encodings 2 columns per field if flag else no
    df_encoded = pd.DataFrame()
    for coln, cold in df.iteritems():
        fe_coln = field_encodings[coln]
        field_encodings_relevant = sorted([int(fe) for fe in fe_coln.keys() if int(fe)>=0])
        # case 1, if there are two types of relevant categories for this field,
        # just return low and high
        if len(field_encodings_relevant) == 2:
            df_encoded[coln] = cold
        # case 2, if there are more than two types of relevant categories for this field
        elif len(field_encodings_relevant) > 2:
            # check if the field is in onesided fields
            if coln in onesided_fields:
                # ignore the first encoding assuming it's answer is no and take all the other encodings into account
                # example bipolar disorder: highs include all types of bipolar disorder
                fe_highs = field_encodings_relevant[1:]
                fe_high_val = "|".join(["-".join(fe_coln[str(fe_high)].replace(",", "").split()) for fe_high in fe_highs])          
                binarized_high = cold.isin(fe_highs).astype(int) 
                df_encoded[f"{coln}_{fe_high_val}_high"] = binarized_high              
            else:
                # check how many low bins will give us at least 25% samples 
                # excluding the highest bin
                for lower_bin_range in range(1, len(field_encodings_relevant)):
                    fe_lows = field_encodings_relevant[:lower_bin_range]
                    fe_low_val = "|".join(["-".join(fe_coln[str(fe_low)].replace(",", "").split()) for fe_low in fe_lows])
                    binarized_low = cold.isin(fe_lows).astype(int)
                    if (sum(binarized_low)/len(binarized_low)) >= 0.25:
                        break        
                # check how many high bins will give us at least 25% samples 
                # exclude the lowest bin
                for higher_bin_range in range(len(field_encodings_relevant) - 1, lower_bin_range - 1, -1):
                    fe_highs = field_encodings_relevant[higher_bin_range:]
                    fe_high_val = "|".join(["-".join(fe_coln[str(fe_high)].replace(",", "").split()) for fe_high in fe_highs])          
                    binarized_high = cold.isin(fe_highs).astype(int)
                    if (sum(binarized_high)/len(binarized_high)) >= 0.25:
                        break
                # check to see if binarized low is the exact opposite of binarized high, then only keep binarized high
                if (binarized_low!=1).astype(int).equals(binarized_high):
                    df_encoded[f"{coln}_{fe_high_val}_high"] = binarized_high
                else:
                    df_encoded[f"{coln}_{fe_low_val}_low"] = binarized_low
                    df_encoded[f"{coln}_{fe_high_val}_high"] = binarized_high
    return df_encoded

def normalize_catsingle(catsingle_df, sel_fields, field_encodings, onesided_fields):
    # select relevant fields
    selected_df = catsingle_df.loc[:, sel_fields]
    # get rid of negative values
    selected_noneg_df = replace_negs(selected_df).dropna(axis=0)
    # get dummies
    selected_noneg_dummies_df = normalize_catsingle_helper(selected_noneg_df, field_encodings, onesided_fields) #pd.get_dummies(selected_noneg_df, columns=sel_fields, drop_first=True)
    # get pca components
    if selected_noneg_dummies_df.shape[1] == 1:
        pca_df = normalize_single_field_function(selected_noneg_dummies_df)
    else:
        pca_df = feature_pca(selected_noneg_dummies_df)    
    return pca_df

def normalize(field_name, sel_fields, field_type_dict, onesided_fields):
    if "integer" in field_name:
        return normalize_integers(field_type_dict["integer"][0], sel_fields)
    
    elif "catsingle" in field_name:
        return normalize_catsingle(field_type_dict["catsingle"][0], sel_fields, field_type_dict["catsingle"][1], onesided_fields)
    
    elif "catmultiple" in field_name:
        return normalize_catsingle(field_type_dict["catmultiple"], sel_fields, onesided_fields)

    return

##########
# encode #
##########

def create_labels(df, lowq=0.25, highq=0.75):
    binary_columns = list(df.columns[df.isin([0, 1, np.nan]).all()])
    non_binary_columns = [c for c in df.columns if c not in binary_columns]
    non_binary_df = df.loc[:, non_binary_columns]
    qdf = non_binary_df.quantile([lowq, highq])
    label_df = non_binary_df.copy()
    label_df = label_df.mask(non_binary_df<qdf.iloc[0, :], other="low") # stringent cutoff for the low field because often these are normal
    label_df = label_df.mask(non_binary_df>=qdf.iloc[1, :], other="high")
    label_df = label_df.mask((non_binary_df>=qdf.iloc[0, :])&(non_binary_df<qdf.iloc[1, :]), other="medium")
    return label_df, df.loc[:, binary_columns]


############
# rarecomb #
############

def run_rarecomb(rare_in, pentities, sentities, rare_out, combos=2):
    cmd = [
        "bash", f"{CURRENT_DIR_PATH}/scripts/run_rarecomb.sh", 
        rare_in, pentities, sentities, rare_out, str(combos)
        ]
    subprocess.run(cmd)
    return
