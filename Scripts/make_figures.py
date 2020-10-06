#!/usr/bin/env python
# coding: utf-8
"""
Authors:
  - Snorre Sulheim, snorres.sulheim@sintef.no
  - Fredrik Fossheim

Date: 17.09.2020
Lisence: CC-BY-4.0

This is the main file used to create figure panels for manuscript.

"""
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import squarify
# plt.style.use("seaborn")
from Bio import SeqIO
import matplotlib.ticker as mtick
from pathlib import Path
import cobra


plt.rcParams.update({'font.size': 22})

prediction_scores = {
    "Domain function \nin total": [177, 215],
    "Extender units": [72, 92],
    "Non-extending \ndomains": [103, 123],
    "Starter units": [5, 8]}
RiPPS = ["lanthipeptide", "thiopeptide", "lassopeptide"]


def figure_prediction_accuracy_bar_chart():
    fig, ax = plt.subplots(1, figsize = (10, 6))
    df = pd.DataFrame(prediction_scores).T
    df.columns = ["Correct", "Total"]
    df["Ratio"] = 100*df["Correct"] / df["Total"]
    print(df)
    b = sns.barplot(y = df.index, x = "Ratio", data = df, color = "#DE4F3C", alpha = 0.8)
    sns.despine()
    ax.set_xlabel("Correct predictions [%]")
    ax.set_ylabel("")

    plt.subplots_adjust(left = 0.4, bottom = 0.2)
    # df["Ratio"].plot(use_index = True, kind = "barh", ax = ax)
    plt.savefig("Figure1c_v1.svg")
    # plt.show()
    print(df)

def prediction_accuracy_stacked_bar_chart():
    

def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

def _rename_BGC_class(df, key = "BGC type"):
    bgc_type_merge =  []
    for x in  df[key]:
        lst = x.split("-")
        if (len(lst) > 1) and lst[0]==lst[1]:
            bgc_class = lst[0]
        else:
            bgc_class = x
        if bgc_class.strip() in RiPPS:
            bgc_class = "RiPPs"
        bgc_type_merge.append(bgc_class)
    return bgc_type_merge

def figure_bigmec_coverage_pie():
    """ Create the figure displaying the number of successful bigmec constructions"""
    filename = "../Data/constructed_pathways/summary.csv"
    df = pd.read_csv(filename, index_col = 0)
    # Rename

    df_success = df.loc[df["Success"]== 1, :]


    df_success["BGC class"] = _rename_BGC_class(df_success)

    cluster_types = df_success.groupby("BGC class").count()
    cluster_types["BGC"].name = ""
    # Convert to treemap
    fig = cluster_types["BGC"].plot(title = None, figsize = (6, 6), kind = "pie")
    plt.savefig("../Figures/bigmec_pie_chart.svg")

def figure_bigmec_coverage_treemap():
    """ Create the figure displaying the number of successful bigmec constructions"""
    filename = "../Data/constructed_pathways/summary.csv"
    df = pd.read_csv(filename, index_col = 0)
    # Rename
    df_success = df.loc[df["Success"]== 1, :]
    bgc_type_merge =  []
    for x in  df_success["BGC type"]:
        lst = x.split("-")
        if (len(lst) > 1) and lst[0]==lst[1]:
            bgc_type_merge.append(lst[0])
        else:
            bgc_type_merge.append(x)
            
    df_success["BGC class"] = bgc_type_merge
    cluster_types = df_success.groupby("BGC class").count()
    cluster_types["BGC"].name = ""
    # Convert to treemap
    labels = []
    for i, row in cluster_types.iterrows():
        labels.append("{0}\n{1}".format(i, row["BGC"]))

    fig, ax = plt.subplots(1, figsize = (10,10))
    l = squarify.plot(label = cluster_types.index, sizes = cluster_types["BGC"], value = cluster_types["BGC"], color=sns.color_palette("tab10"), ax = ax)
    plt.axis("off")
    plt.savefig("../Figures/bigmec_treemap.svg")


def figure_prediction_coverage():
    fn = "../Data/knockouts/all_optknock_results_0.5.csv"
    df = pd.read_csv(fn)

    print(df)

def figure_prediction_reactions():
    fn = "../Data/knockouts/all_optknock_results_0.5.csv"
    df = pd.read_csv(fn)
    
    count_df = df.groupby(['ID', 'BGC type']).size().unstack(fill_value=0)
    count_df_ratio = count_df/count_df.sum(axis = 0)
    count_df_ratio.T.plot(kind = "barh", stacked = True, cmap = "tab20", figsize = (12,6))
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))    
    plt.subplots_adjust(left = 0.3, right = 0.7)
    plt.savefig("../Figures/prediction_reactions.svg")

def figure_prediction_reactions_max():
    fn = "../Data/knockouts/all_optknock_results_0.5.csv"
    df = pd.read_csv(fn)
    g = df.groupby("BGC")
    max_rows = []
    for i, group in g:
        maxx = group["Production"].max()
        ix = group["Production"] == maxx
        selected = group.loc[ix, :]
        max_rows += list(selected.index)
    print(max_rows)
    max_df = df.loc[max_rows, ["ID", "BGC type"]]
    count_df = max_df.groupby(['ID', 'BGC type']).size().unstack(fill_value=0)
    count_df_ratio = 100*count_df/count_df.sum(axis = 0)
    ax = count_df_ratio.T.plot(kind = "barh", stacked = True, cmap = "tab20", figsize = (12, 6))
    ax.xaxis.set_major_formatter(mtick.PercentFormatter())
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))    
    plt.subplots_adjust(left = 0.3, right = 0.7)
    plt.savefig("../Figures/prediction_reactions_max.svg")


def plot_taxonomic_pie():
    filename = "../Data/constructed_pathways/summary.csv"
    df = pd.read_csv(filename, index_col = 0)
    # Rename
    df_success = df.loc[df["Success"]== 1, :]

    bgc_ids = list(df_success["BGC"])
    taxonomy_list = []
    for i in bgc_ids:
        fn = "../Data/mibig/{0}.gbk".format(i)
        tax = _get_strain_knowledge(fn)
        taxonomy_list.append(tax)

    df_tax = pd.DataFrame(taxonomy_list)
    df_tax = df_tax.sort_values(0)
    g1 = df_tax.groupby(1).size()
    #g01 = df_tax.groupby([0, 1]).size().unstack(fill_value = 0)
    g1["Others"] = g1[g1<5].sum()
    g1.name = ""    
    fig, ax = plt.subplots(1, figsize = (12, 12))
    plt.set_cmap("tab10")
    pie,_,txt = ax.pie(g1[g1>5], radius = 1,autopct='%1.2f%%',  pctdistance=.8)
    width = 0.4
    plt.setp(pie, width=width, edgecolor='white')
    plt.legend(labels = g1[g1>5].index, loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.subplots_adjust(left = 0.3, right = 0.7)
    plt.savefig("../Figures/taxonomy_pie.svg")

def plot_taxonomic_treemap():
    filename = "../Data/constructed_pathways/summary.csv"
    df = pd.read_csv(filename, index_col = 0)
    # Rename
    df_success = df.loc[df["Success"]== 1, :]

    bgc_ids = list(df_success["BGC"])
    taxonomy_list = []
    for i in bgc_ids:
        fn = "../Data/mibig/{0}.gbk".format(i)
        tax = _get_strain_knowledge(fn)
        taxonomy_list.append(tax)

    df_tax = pd.DataFrame(taxonomy_list)
    df_tax = df_tax.sort_values(0)
    g1 = df_tax.groupby(1).size()
    #g01 = df_tax.groupby([0, 1]).size().unstack(fill_value = 0)
    g1["Others"] = g1[g1<5].sum()
    g1.name = ""    
    fig, ax = plt.subplots(1, figsize = (10, 10))
    #pie,_,txt = ax.pie(g1[g1>5], radius = 1,autopct='%1.2f%%',  pctdistance=.8)
    squarify.plot(label = g1[g1>5].index, sizes = g1[g1>5], value = g1[g1>5], norm_x = 30, color = sns.color_palette("tab10"), ax = ax)
    plt.axis("off")
    plt.savefig("../Figures/taxonomy_treemap.svg")
    #plt.setp(txt, size = 12)#, weight="bold")

def _get_strain_knowledge(cluster_path):  # yes
    gb_list = []
    for gb_record in SeqIO.parse(open(cluster_path, "r"), "genbank"):
        gb_list.append(gb_record)
    g = gb_list[0]
    return g.annotations["taxonomy"]


def pathway_length_figure():
    folder = Path("../Data/constructed_pathways/")
    model_info_list = []
    substrate_dct_list = []
    for fn in folder.glob("*.json"):
        model = cobra.io.load_json_model(str(fn))
        n_reactions_total, n_reactions_link, substrate_dct = _get_model_info(model)
        model_info_list.append([fn.stem, n_reactions_link, n_reactions_total])
        substrate_dct_list.append(substrate_dct)

    df_info = pd.DataFrame(model_info_list, columns = ["ID", "Pathway length", "Total number of reactions"])
    df_info["ID"] = df_info["ID"].astype(int)

    df_success = _get_df_success()
    print(df_info)
    print(df_success)
    df_info = df_info.merge(df_success[["BGC", "BGC class"]], left_on="ID", right_on= "BGC")
    sns.histplot(data = df_info[df_info["Pathway length"]>1], x = "Pathway length", hue = "BGC class", 
             kde = True, element = "step")
    plt.savefig("../Figures/pathway_length.svg")

    df_substrates = df_substrates.abs()

def _get_df_success():
    filename = "../Data/constructed_pathways/summary.csv"
    df = pd.read_csv(filename, index_col = 0)
    # Rename
    df_success = df.loc[df["Success"]== 1, :]
    return df_success

def _get_model_info(model):
    #    for r in model.reactions:
    #        print(r.id, "\t", r.reaction)
    n_reactions_total = len(model.reactions)
    n_reactions_link = 0
    for r in model.reactions:
        
        try:
            i = int(r.id.split("_")[-1])
        except:
            pass
        else:
            if i > n_reactions_link:
                n_reactions_link = i
    
    lump_reaction = _lump_reaction(model)
    substrate_dct = _get_substrate_dict(lump_reaction)
    return n_reactions_total, n_reactions_link, substrate_dct
    
def _lump_reaction(model):
    for i, r in enumerate(model.reactions):
        if i == 0:
            lump_reaction = r
        else:
            lump_reaction += r
    return lump_reaction

def _get_substrate_dict(r):
    dct = {}
    for m, i in r.metabolites.items():
        if i < 0:
            dct[m.id]=i
    return dct
    

if __name__ == '__main__':
    if 0:
        figure_prediction_accuracy_bar_chart()
    if 0:
        figure_bigmec_coverage_pie()
        figure_bigmec_coverage_treemap()
    if 0:
        figure_prediction_coverage()

    if 0:
        figure_prediction_reactions()
        figure_prediction_reactions_max()

    if 0:
        # plot_taxonomic_pie()
        plot_taxonomic_treemap()

    if 1:
        pathway_length_figure()
