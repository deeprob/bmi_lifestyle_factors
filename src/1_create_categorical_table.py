import pandas as pd
import numpy as np
import os
import json
import utils as ut



def main(lifestyle_df, input_dir, field_encoding_dir, json_files, onesided_fields_files, store_dirs):
    int_df, int_fe_dicts = ut.get_integer_factors(lifestyle_df, input_dir, field_encoding_dir)
    catsingle_df, catsingle_fe_dicts = ut.get_catsingle_factors(lifestyle_df, input_dir, field_encoding_dir)
    catmultiple_df, catmultiple_fe_dicts = ut.get_catmultiple_factors(lifestyle_df, input_dir, field_encoding_dir)

    field_type_dict = {
        "integer": (int_df, int_fe_dicts),
        "catsingle": (catsingle_df, catsingle_fe_dicts),
        "catmultiple": (catmultiple_df, catmultiple_fe_dicts)
    }
    for json_file, onesided_fields_file, store_dir in zip(json_files, onesided_fields_files, store_dirs):
        with open(json_file, "r") as f:
            cat_dict = json.load(f)
        
        onesided_fields = []
        with open(onesided_fields_file) as f:
            for lines in f:
                onesided_fields.append(int(lines.strip()))

        # normalize and decompose fields
        pca_raw_dfs = [ut.normalize(fn, sfields, field_type_dict, onesided_fields) for fn,sfields in cat_dict.items()]
        # save tables
        for fn, pca_raw_df in zip(cat_dict.keys(), pca_raw_dfs):
            pca_raw_df.to_csv(os.path.join(store_dir, f"{fn}.csv"))

        # create cat table
        cat_table = pd.concat([df.loc[:, "component_1"].rename(coln) for df,coln in zip(pca_raw_dfs, cat_dict.keys())], axis=1)
        cat_table.to_csv(os.path.join(store_dir, "cat.csv"))
        # label cat table excluding the fields which are already boolean encoded
        cat_table_labelled, binary_table = ut.create_labels(cat_table)
        cat_table_labelled.to_csv(os.path.join(store_dir, "cat_labelled.csv"))
        # one hot encode cat table
        cat_table_ohe = pd.get_dummies(cat_table_labelled, dummy_na=True)
        # add high to binary table
        binary_table.columns = [f"{c}_high" for c in binary_table.columns]
        # fill binary table nan values with 0
        binary_table = binary_table.fillna(0)
        cat_table_ohe = pd.concat((cat_table_ohe, binary_table), axis=1)
        cat_table_ohe.to_csv(os.path.join(store_dir, "cat_encoded.csv"))

        # keep only high low :TODO: filter fields with less than 1% samples with a 1 encoding 
        cat_table_selected = cat_table_ohe.loc[:, [c for c in cat_table_ohe.columns if (c.endswith("_high") or c.endswith("_low"))]]
        cat_table_selected.columns = [f"Input_{c}" for c in cat_table_selected.columns]
        cat_table_selected.to_csv(os.path.join(store_dir, "cat_extremes.csv"))
    return

if __name__ == "__main__":
    lifestyle_file = "/data5/deepro/ukbiobank/analysis/lifestyle_factors/data/lifestyle_v2.xlsx"
    json_files = [
        "/data5/deepro/ukbiobank/analysis/lifestyle_factors/data/mental_health/fields.json",
        "/data5/deepro/ukbiobank/analysis/lifestyle_factors/data/physical_activity/fields.json",
        "/data5/deepro/ukbiobank/analysis/lifestyle_factors/data/diet/fields.json",
        ]
    merged_field_input_dir = "/data5/deepro/ukbiobank/preprocess/rarecomb_pheno_prepare/data"
    field_encoding_dir = "/data5/deepro/ukbiobank/download/download_phenotypes/data"
    store_dirs = [
        "/data5/deepro/ukbiobank/analysis/lifestyle_factors/data/mental_health",
        "/data5/deepro/ukbiobank/analysis/lifestyle_factors/data/physical_activity",
        "/data5/deepro/ukbiobank/analysis/lifestyle_factors/data/diet",
        ]
    onesided_fields_files = [
        "/data5/deepro/ukbiobank/analysis/lifestyle_factors/data/mental_health/onesided.txt",
        "/data5/deepro/ukbiobank/analysis/lifestyle_factors/data/physical_activity/onesided.txt",
        "/data5/deepro/ukbiobank/analysis/lifestyle_factors/data/diet/onesided.txt",
        ]

    lifestyle_df = pd.read_excel(lifestyle_file)
    # take one sided file as input
    main(lifestyle_df, merged_field_input_dir, field_encoding_dir, json_files, onesided_fields_files, store_dirs)

