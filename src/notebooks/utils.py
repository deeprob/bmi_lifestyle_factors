import os
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.pyplot as plt
import prince

sns.set_style("dark", {"grid.color": ".6", "grid.linestyle": ":"})
sns.set_context("paper")
sns.set(font_scale = 1.5)

##############################
# read the lifestyle factors #
##############################

def get_continuous_factors(lifestyle_df, store_dir):
    lifestyle_df = lifestyle_df.loc[(lifestyle_df.shortlist=="X") & (lifestyle_df.Type=="continuous")]
    dfs = [pd.read_csv(os.path.join(store_dir, t, g, "tables", f"{i}_quantile.csv"), usecols=["eid", "merged"], index_col=["eid"]).rename(columns={"merged": i}) for t,g,i in zip(lifestyle_df.Type, lifestyle_df.Phenotype_group, lifestyle_df.Phenotype_ID)]
    df = pd.concat(dfs, axis=1).sort_values("eid")
    return df


def get_integer_factors(lifestyle_df, store_dir):
    lifestyle_df = lifestyle_df.loc[(lifestyle_df.shortlist=="X") & (lifestyle_df.Type=="integer")]
    dfs = [pd.read_csv(os.path.join(store_dir, t, g, "tables", f"{i}_quantile.csv"), usecols=["eid", "merged"], index_col=["eid"]).rename(columns={"merged": i}) for t,g,i in zip(lifestyle_df.Type, lifestyle_df.Phenotype_group, lifestyle_df.Phenotype_ID)]
    df = pd.concat(dfs, axis=1).sort_values("eid")
    return df

def get_catsingle_factors(lifestyle_df, store_dir):
    lifestyle_df = lifestyle_df.loc[(lifestyle_df.shortlist=="X") & (lifestyle_df.Type=="categorical_single")]
    dfs = [pd.read_csv(os.path.join(store_dir, t, g, "tables", f"{i}.csv"), usecols=["eid", "merged"], index_col=["eid"]).rename(columns={"merged": i}) for t,g,i in zip(lifestyle_df.Type, lifestyle_df.Phenotype_group, lifestyle_df.Phenotype_ID)]
    df = pd.concat(dfs, axis=1).sort_values("eid")
    return df

def get_catmultiple_factors(lifestyle_df, store_dir):
    lifestyle_df = lifestyle_df.loc[(lifestyle_df.shortlist=="X") & (lifestyle_df.Type=="categorical_multiple")]
    dfs = [pd.read_csv(os.path.join(store_dir, t, g, "tables", f"{i}.csv"), usecols=lambda x: x.startswith("eid") or x.startswith("binarized"), index_col=["eid"]).rename(columns=lambda col: col.replace("binarized", f"{i}")) for t,g,i in zip(lifestyle_df.Type, lifestyle_df.Phenotype_group, lifestyle_df.Phenotype_ID)]
    df = pd.concat(dfs, axis=1).sort_values("eid")
    return df

#############
# normalize #
#############

def standardize_df(df):
    # df = (df - df.mean())/df.std()
    df = (df-df.min())/(df.max()-df.min())
    return df

def replace_negs(df):
    return df.mask(df<0)
    
def create_labels(df, lowq=0.25, highq=0.75):
    binary_columns = list(df.columns[df.isin([0, 1, np.nan]).all()])
    non_binary_columns = [c for c in df.columns if c not in binary_columns]
    non_binary_df = df.loc[:, nonbinary_columns]
    qdf = non_binary_df.quantile([lowq, highq])
    label_df = non_binary_df.copy()
    label_df = label_df.mask(non_binary_df<qdf.iloc[0, :], other="low")
    label_df = label_df.mask(non_binary_df>qdf.iloc[1, :], other="high")
    label_df = label_df.mask((non_binary_df>=qdf.iloc[0, :])&(non_binary_df<=qdf.iloc[1, :]), other="medium")
    return label_df, df.loc[:, binary_columns]

def visualize_field_pca(df):
    pca = PCA(n_components=2)
    Xt = pca.fit_transform(df.values)
    plot_df = pd.DataFrame()
    plot_df["fields"] = df.index
    plot_df["component_1"] = Xt[:,0]
    plot_df["component_2"] = Xt[:,1]
    fig, axes = plt.subplots(figsize=(12,10))
    sns_ax = sns.scatterplot(
        data=plot_df, 
        x="component_1", y="component_2", 
        hue="fields", style="fields", 
        s=200, palette="Set1", 
        ax=axes)
    return sns_ax

def feature_pca(df, label_df):
    pca = PCA(n_components=2)
    Xt = pca.fit_transform(df.values)
    plot_df = pd.DataFrame()
    plot_df.index = df.index
    plot_df["component_1"] = Xt[:,0]
    plot_df["component_2"] = Xt[:,1]
    plot_df = plot_df.merge(label_df, left_index=True, right_index=True)
    return plot_df

def feature_mca(df):
    mca = prince.MCA(
        n_components=2, 
        n_iter=3, 
        copy=True, 
        check_input=True, 
        engine='auto', 
        random_state=42
        )
    print(df.loc[np.isinf(df).any(axis=1)])
    print(df.loc[pd.isnull(df).any(axis=1)])
    mca = mca.fit(df)
    plot_df = mca.row_coordinates(df)
    plot_df.index = df.index
    plot_df.columns = ["component_1", "component_2"]
    plot_df = plot_df.merge(df, left_index=True, right_index=True)
    return plot_df
