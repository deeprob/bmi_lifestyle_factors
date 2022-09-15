import os
import pandas as pd


def main(in_dir, in_sub_dir, infilename, wes_table, cases, controls):

    # read wes file
    wes_df = pd.read_csv(wes_table, sep="\t", index_col=0)
    print("wes loaded")
    # read case eids
    cases_list = []
    with open(cases, "r") as f:
        for lines in f:
            cases_list.append(int(lines.strip()))
    # read controls eids
    controls_list = []
    with open(controls, "r") as f:
        for lines in f:
            controls_list.append(int(lines.strip()))

    # create categorical tables for rarecomb input
    lifestyle_dfs = []
    for isd in in_sub_dir:
        sub_dir = os.path.join(in_dir, isd)
        infilepath = os.path.join(sub_dir, infilename)
        # read encoded extreme lifestyle factors
        lifestyle_df = pd.read_csv(infilepath, index_col=0)
        lifestyle_dfs.append(lifestyle_df)
        # merge lifestyle factors with wes data
        df = wes_df.merge(lifestyle_df, left_index=True, right_index=True)
        df.index.name = "Sample_Name"
        # only keep case and control samples
        df = df.loc[[i for i in cases_list + controls_list if i in df.index]]
        # create outpur pheno column
        df["Output_pheno"] = 0
        df.loc[[i for i in cases_list if i in df.index], "Output_pheno"] = 1 
        # save rarecomb table
        rarecomb_infile = os.path.join(sub_dir, "cat_rarecomb_in.csv")
        df.to_csv(rarecomb_infile)
        # save primary entities
        primary_entities_file = os.path.join(sub_dir, "primary.csv")
        with open(primary_entities_file, "w") as f:
            for c in lifestyle_df.columns:
                f.write(f"{c}\n")
        # save secondary entities
        secondary_entities_file = os.path.join(sub_dir, "secondary.csv")
        with open(secondary_entities_file, "w") as f:
            for c in wes_df.columns:
                f.write(f"{c}\n") 

    # create meta table for rarecomb input
    sub_dir = in_dir
    # concatenate meta lifestyle factors
    lifestyle_df = pd.concat(lifestyle_dfs, axis=1)
    # merge lifestyle factors with wes data
    df = wes_df.merge(lifestyle_df, left_index=True, right_index=True)
    df.index.name = "Sample_Name"
    # only keep case and control samples
    df = df.loc[[i for i in cases_list + controls_list if i in df.index]]
    # fill na values which is a result of pd concat, rarecomb can't work with nas
    df = df.fillna(0)
    # create output pheno column
    df["Output_pheno"] = 0
    df.loc[[i for i in cases_list if i in df.index], "Output_pheno"] = 1 
    # save rarecomb table
    rarecomb_infile = os.path.join(sub_dir, "meta_rarecomb_in.csv")
    df.to_csv(rarecomb_infile)
    # save primary entities
    primary_entities_file = os.path.join(sub_dir, "primary.csv")
    with open(primary_entities_file, "w") as f:
        for c in lifestyle_df.columns:
            f.write(f"{c}\n")
    # save secondary entities
    secondary_entities_file = os.path.join(sub_dir, "secondary.csv")
    with open(secondary_entities_file, "w") as f:
        for c in wes_df.columns:
            f.write(f"{c}\n")           
    return


if __name__ == "__main__":
    wes_table = "/data5/UK_Biobank/bmi_project/combinations/white_british/tables/wes.tsv"
    in_dir = "/data5/deepro/ukbiobank/analysis/lifestyle_factors/data/"
    in_sub_dirs = ["mental_health", "physical_activity", "diet"]
    infilename = "cat_extremes.csv"
    pheno_cases = "/data5/UK_Biobank/bmi_project/combinations/white_british/cases_controls/high_low_cases.txt"
    pheno_controls = "/data5/UK_Biobank/bmi_project/combinations/white_british/cases_controls/high_low_controls.txt"
    main(in_dir, in_sub_dirs, infilename, wes_table, pheno_cases, pheno_controls)
