#!/usr/bin/env python
# coding: utf-8
"""
Copyright 2020 Snorre Sulheim (snorre.sulheim@sintef.no)
https://github.com/AlmaasLab/BiGMeC

This is the file used to make figure panels for the reseach paper associated with BiGMeC.


This file is part of BiGMeC. BiGMeC is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. BiGMeC is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with BiGMeC.
If not, see <http://www.gnu.org/licenses/>.

Authors:
  - Snorre Sulheim, snorres.sulheim@sintef.no
  - Fredrik A. Fossheim

Date: 17.09.2020
"""
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import squarify
plt.style.use("seaborn-white")
from Bio import SeqIO
import matplotlib.ticker as mtick
from pathlib import Path
import cobra
import numpy as np
from scipy.stats import pearsonr

import venn
import catheat

from matplotlib.colors import ListedColormap

# cmap = ListedColormap(sns.color_palette())
# plt.set_cmap(cmap)
# Set the default color cycle
# plt.rcParams['axes.prop_cycle'] = plt.cycler(color=sns.color_palette("muted")) 
sns.set_style("white")
sns.color_palette("muted")

plt.rcParams.update({'font.size': 22})

prediction_scores = {
    "Domain function \nin total": [177, 215],
    "Extender units": [72, 92],
    "Non-extending \ndomains": [105, 123],
    "Starter units": [5, 8]}
RiPPS = ["lanthipeptide", "thiopeptide", "lassopeptide"]

pathway_mibig_dict = {
    "Anabaenopeptin": 302,
    "Bafilomycin": 28,
    "Difficidin": 176,
    "Geldanamycin": 66,
    "Leupyrrin": 380,
    "Oocydin": 1032,
    "Oxazolomycin": 1106,
    "Tolaasin": 447
}

def figure_production_in_real_vs_constructed():
    model_fn = "../Models/Sco-GEM.xml"
    model = cobra.io.read_sbml_model(model_fn)
    print(model.medium)

    folders = ["../Data/validation_pathways/", "../Data/constructed_pathways/"]
    wt_growth = model.slim_optimize()
    model.reactions.BIOMASS_SCO_tRNA.lower_bound = 0.9 * wt_growth
    production_dict = {}
    for name, bgc_int in pathway_mibig_dict.items():
        print(name)
        production_list = []
        mets_list = []
        for i, model_fn in enumerate([name, bgc_int]):
            pathway_fn = folders[i] +  "{0}.json".format(model_fn)
            pathway = cobra.io.load_json_model(pathway_fn)
            with model:
                model.merge(pathway)
                model.objective = model.reactions.DM_secondary_metabolite
                production = model.slim_optimize()
            production_list.append(production)
        production_dict[name] = production_list
    df = pd.DataFrame(production_dict)
    df.index = ["Real", "Constructed from BGC"]
    fig, ax = plt.subplots(1, figsize = (10, 6))
    df.T.plot(kind = "bar", alpha = 0.8, ax = ax)
    ax.set_ylabel("Production rate [mmol/gDW/h]")
    ax.set_xlabel("")
    plt.subplots_adjust(bottom = 0.3)
    plt.savefig("../Figures/production_rate_comparison.svg")

def figure_production_in_real_vs_constructed_scatter():
    model_fn = "../Models/Sco-GEM.xml"
    model = cobra.io.read_sbml_model(model_fn)
    print(model.medium)
    folders = ["../Data/validation_pathways/", "../Data/constructed_pathways/"]
    wt_growth = model.slim_optimize()
    print(model.summary())
    model.reactions.BIOMASS_SCO_tRNA.lower_bound = 0.9 * wt_growth
    production_dict = {}
    for name, bgc_int in pathway_mibig_dict.items():
        print(name)
        production_list = []
        mets_list = []
        for i, model_fn in enumerate([name, bgc_int]):
            pathway_fn = folders[i] +  "{0}.json".format(model_fn)
            pathway = cobra.io.load_json_model(pathway_fn)
            with model:
                model.merge(pathway)
                model.objective = model.reactions.DM_secondary_metabolite
                production = model.slim_optimize()
            production_list.append(production)
        production_dict[name] = production_list
    df = pd.DataFrame(production_dict)
    df.index = ["Real", "Constructed from BGC"]
    fig, ax = plt.subplots(1, figsize = (10, 6))
    ax.plot([0,0.1], [0,0.1], ls = "--", c = "k", alpha = 0.5)
    for i, bgc_name in enumerate(df.columns):
        x = df.loc["Real", bgc_name]
        y = df.loc["Constructed from BGC", bgc_name]
        ax.scatter(x, y, color = "C{0}".format(i), label = bgc_name, alpha = 0.8, s = 200)
    ax.set_ylabel("BiGMeC-constructed \n production rate [mmol/gDW/h]")
    ax.set_xlabel("Literature-constructed \n production rate [mmol/gDW/h]")
    # sns.despine()
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(right = 0.6, bottom = 0.2, left = 0.2)
    ax.set_xlim([1e-3, df.max().max()*1.1])
    ax.set_ylim([1e-3, df.max().max()*1.1])
    plt.savefig("../Figures/production_rate_comparison_scatter.svg")

    # Calculate correlation
    x = df.loc["Real", :].values
    y = df.loc["Constructed from BGC", :].values
    (corr, p_value) = pearsonr(x, y)
    print("Correlation: ", corr)
    print("P-value: ", p_value)

def figure_knockouts_predictions_in_real_vs_constructed():
    plt.rcParams.update({'font.size': 18})
    folders = ["../Data/validation_pathways/knockouts", "../Data/knockouts"]
    knockout_dict = {}
    for name, bgc_int in pathway_mibig_dict.items():
        real_df = pd.read_csv(folders[0]+"/brute_optknock_{0}_0.5.csv".format(name))
        constructed_df = pd.read_csv(folders[1]+"/brute_optknock_{0}_0.5.csv".format(bgc_int))
        reactions = list(real_df["ID"])+list(constructed_df["ID"])
        if len(real_df["ID"])==0:
            reactions.append("None")
        if len(constructed_df["ID"])==0:
            reactions.append("None")
        r_dict = {}
        for r in set(reactions):
            r_dict[r] = str(int(reactions.count(r)))
        knockout_dict[name] = r_dict

    df = pd.DataFrame(knockout_dict).T
    df["None"] = df.pop("None")
    mask = df.isna()
    df = df.astype(str)
    fig, ax = plt.subplots(1, figsize = (10, 6))
    ax = catheat.heatmap(df, square = True, linewidths=2, linecolor='w', cmap = {"nan":"#e5e5e5", "1":"C3","2":"C2"}, palette = "muted", alpha = 0.8, xticklabels=True, yticklabels=True, ax = ax)#, mask = mask)
    #ax[0].set_facecolor('#e5e5e5')
    plt.subplots_adjust(left = 0.3, right = 0.9, bottom = 0.25) 
    plt.savefig("../Figures/production_rate_comparison_knockout.svg")           

def figure_prediction_accuracy_bar_chart():
    fig, ax = plt.subplots(1, figsize = (10, 6))
    df = pd.DataFrame(prediction_scores).T
    df.columns = ["Correct", "Total"]
    df["Ratio"] = 100*df["Correct"] / df["Total"]

    print(df)
    b = sns.barplot(y = df.index, x = "Ratio", data = df, alpha = 0.8)
    sns.despine()
    ax.set_xlabel("Correct predictions [%]")
    ax.set_ylabel("")
    ax.set_xlim(0, 100)
    autolabel(b.patches, ax, df)
    plt.subplots_adjust(left = 0.4, bottom = 0.2)
    # df["Ratio"].plot(use_index = True, kind = "barh", ax = ax)
    plt.savefig("../Figures/prediction_accuracy_full.svg")
    # plt.show()
    #print(df)

def autolabel(rects, ax, df):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for i, rect in enumerate(rects):
        print(rect.get_height(), rect.get_y(), rect.get_y()-rect.get_height()/2)
        ax.annotate('{0}/{1}'.format(df["Correct"][i], df["Total"][i]),
                    xy=(rect.get_width(), rect.get_y()+rect.get_height()/2),    
                    xytext=(2, 0),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='left', va='center')

def figure_prediction_accuracy_full_vertical():
    fn = "../Data/prediction_accuracy_8_experimental.xlsx"
    df = pd.read_excel(fn, index_col = 0)
    df.drop(columns = "Note", index = "SUM", inplace = True)
    df.sort_index(inplace = True, ascending = True)
    fig, ax = plt.subplots(1, figsize = (16, 8))
    width = 0.2
    space = 0.26
    x = np.arange(8)-(width+space)/2

    bars = ax.bar(x, df.loc[:, "Total domains"], width, color = "none", edgecolor = "C0", alpha = 0.8, label = "Total domains")
    bars2 = ax.bar(x, df.loc[:, "Correct domains"],width,  edgecolor = "C0", alpha = 0.8)
    for bar in bars:
        bar.set_edgecolor("C0")
        bar.set_linewidth(2)
    autolabel2_vertical(bars, bars2, ax)

    x = np.arange(8)
    bars = ax.bar(x, df.loc[:, "Total KS domains"], width, color = "none", edgecolor = "C1", alpha = 0.8, label = "Extending domains")
    bars2 = ax.bar(x, df.loc[:, "Correct KS domains"], width, edgecolor = "C1", alpha = 0.8)
    for bar in bars:
        bar.set_edgecolor("C1")
        bar.set_linewidth(2)
        
    autolabel2_vertical(bars, bars2, ax)


    x = np.arange(8)+(width+space)/2
    bars = ax.bar(x, df.loc[:, "Total other domains"], width, color = "none", edgecolor = "C2", alpha = 0.8, label = "Non-extending domains")
    bars2 = ax.bar(x, df.loc[:, "Correct other domains"], width, edgecolor = "C2", alpha = 0.8)
    for bar in bars:
        bar.set_edgecolor("C2")
        bar.set_linewidth(2)
        
    autolabel2_vertical(bars, bars2, ax)
    ax.set_xticks(np.arange(len(df.index)+1))
    ax.set_xticklabels(list(df.index), rotation = 30)
    ax.set_ylim(0, 60)
    print(df.index)
    ax.set_ylabel("Number of domains")
    sns.despine()
    plt.subplots_adjust(bottom = 0.3)
    plt.legend()
    plt.savefig("../Figures/bigmec_accuracy_full_vertical.svg")


def autolabel2_vertical(rects, rects2, ax):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect, rect2 in zip(rects, rects2):
        height = rect.get_height()
        ax.annotate('{0}/{1}'.format(rect2.get_height(), height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(1, 5),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', rotation = 90, fontsize = 14)

def figure_prediction_accuracy_full_horizontal():
    fn = "../Data/prediction_accuracy_8_experimental.xlsx"
    df = pd.read_excel(fn, index_col = 0)
    df.drop(columns = "Note", index = "SUM", inplace = True)
    df.sort_index(inplace = True, ascending = False)
    fig, ax = plt.subplots(1, figsize = (10, 8))
    width = 0.2
    space = 0.26
    x = np.arange(8)+(width+space)/2

    bars = ax.barh(x, df.loc[:, "Total domains"], width, color = "none", edgecolor = "C0", alpha = 0.8, label = "Total domains")
    bars2 = ax.barh(x, df.loc[:, "Correct domains"],width,  edgecolor = "C0", alpha = 0.8)
    for bar in bars:
        bar.set_edgecolor("C0")
        bar.set_linewidth(2)
    autolabel2_horizontal(bars, bars2, ax)

    x = np.arange(8)
    bars = ax.barh(x, df.loc[:, "Total KS domains"], width, color = "none", edgecolor = "C1", alpha = 0.8, label = "Extending domains")
    bars2 = ax.barh(x, df.loc[:, "Correct KS domains"], width, edgecolor = "C1", alpha = 0.8)
    for bar in bars:
        bar.set_edgecolor("C1")
        bar.set_linewidth(2)
        
    autolabel2_horizontal(bars, bars2, ax)


    x = np.arange(8)-(width+space)/2
    bars = ax.barh(x, df.loc[:, "Total other domains"], width, color = "none", edgecolor = "C2", alpha = 0.8, label = "Non-extending domains")
    bars2 = ax.barh(x, df.loc[:, "Correct other domains"], width, edgecolor = "C2", alpha = 0.8)
    for bar in bars:
        bar.set_edgecolor("C2")
        bar.set_linewidth(2)
        
    autolabel2_horizontal(bars, bars2, ax)
    ax.set_yticklabels([1]+list(df.index))#, rotation = 90)
    print(df.index)
    ax.set_xlim(0, 52)
    ax.set_xlabel("Number of domains")
    sns.despine()
    plt.legend()
    plt.subplots_adjust(left = 0.4)
    plt.savefig("../Figures/bigmec_accuracy_full_horizontal.svg")

def autolabel2_horizontal(rects, rects2, ax):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect, rect2 in zip(rects, rects2):
        height = rect.get_width()
        ax.annotate('{0}/{1}'.format(rect2.get_width(), height),
                    xy=(rect.get_width(), rect.get_y()+rect.get_height()/2),
                    xytext=(4,-1),#(0, 5),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='left', va='center', fontsize = 11)#, rotation = 90)

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


def figure_bigmec_coverage_venn():
    """ Create the figure displaying the number of successful bigmec constructions"""
    filename = "../Data/constructed_pathways/summary.csv"
    df = pd.read_csv(filename, index_col = 0)
    df_success = df.loc[df["Success"]== 1, :]
    print(df_success.shape)

    dict_list = []
    for i, row in df_success.iterrows():
        dic = {}
        if isinstance(row["BGC type"], float):
            continue
        bgc_types = row["BGC type"].split("/")
        bgc_types
        for k in bgc_types:
            dic[k] = True
        dic["BGC"] = row["BGC"]
        dict_list.append(dic)

    df_venn = pd.DataFrame(dict_list)
    df_venn = df_venn.fillna(False)
    cols = ["T1PKS", "transAT-PKS-like", "transAT-PKS", "NRPS-like", "NRPS", "PKS-like", "BGC"]
    other_columns = [x for x in df_venn.columns if not x in cols]
    df_venn2 = pd.DataFrame()
    df_venn2["Type 1 PKS"] = df_venn["T1PKS"]
    df_venn2["trans-AT PKS"] = df_venn[["transAT-PKS-like", "transAT-PKS"]].sum(axis = 1).astype(bool)
    df_venn2["NRPS"] = df_venn[["NRPS-like", "NRPS"]].sum(axis = 1).astype(bool)
    df_venn2["Other"] = df_venn[other_columns].sum(axis = 1).astype(bool)
    df_venn2.index = df_venn["BGC"]

    dic = {"Type 1 PKS": set(df_venn2[df_venn2["Type 1 PKS"]].index),
           "trans-AT PKS": set(df_venn2[df_venn2["trans-AT PKS"]].index),
           "NRPS": set(df_venn2[df_venn2["NRPS"]].index),
           "Other": set(df_venn2[df_venn2["Other"]].index)}
    cmap = ListedColormap(sns.color_palette("muted"), N = 4)
    print(cmap)
    venn.venn(dic, cmap = cmap, fontsize = 26)
    plt.savefig("../Figures/bigmec_venn.svg")

    print("N total successful: ", len(df_success))
    for key, value in dic.items():
        print(key, len(value))


def figure_bigmec_unsuccessful_coverage_venn():
    """ Create the figure displaying the number of successful bigmec constructions"""
    filename = "../Data/constructed_pathways/summary.csv"
    df = pd.read_csv(filename, index_col = 0)
    df_success = df.loc[df["Success"]== 0, :]

    dict_list = []
    for i, row in df_success.iterrows():
        dic = {}
        if isinstance(row["BGC type"], float):
            continue
        bgc_types = row["BGC type"].split("/")
        bgc_types
        for k in bgc_types:
            dic[k] = True
        dic["BGC"] = row["BGC"]
        dict_list.append(dic)

    df_venn = pd.DataFrame(dict_list)
    df_venn = df_venn.fillna(False)

    print(df_venn.sum())

    cols = ["T1PKS", "transAT-PKS-like", "transAT-PKS", "NRPS-like", "NRPS", "PKS-like", "BGC"]
    other_columns = [x for x in df_venn.columns if not x in cols]
    df_venn2 = pd.DataFrame()
    df_venn2["T1PKS"] = df_venn["T1PKS"]
    df_venn2["TransAT-PKS"] = df_venn[["transAT-PKS-like", "transAT-PKS"]].sum(axis = 1).astype(bool)
    df_venn2["NRPS"] = df_venn[["NRPS-like", "NRPS"]].sum(axis = 1).astype(bool)
    df_venn2["Other"] = df_venn[other_columns].sum(axis = 1).astype(bool)
    df_venn2.index = df_venn["BGC"]

    dic = {"T1PKS": set(df_venn2[df_venn2["T1PKS"]].index),
           "TransAT-PKS": set(df_venn2[df_venn2["TransAT-PKS"]].index),
           "NRPS": set(df_venn2[df_venn2["NRPS"]].index),
           "Other": set(df_venn2[df_venn2["Other"]].index)}
    venn.venn(dic)
    plt.savefig("../Figures/bigmec_venn_unsuccessful.svg")

    print("N total successful: ", len(df_success))
    for key, value in dic.items():
        print(key, len(value))

def figure_bigmec_coverage_treemap():
    """ Create the figure displaying the number of successful bigmec constructions"""
    filename = "../Data/constructed_pathways/summary.csv"
    df = pd.read_csv(filename, index_col = 0)
    # Rename
    df_success = df.loc[df["Success"]== 1, :]
    print("##")
    print(df_success.shape)
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
    l = squarify.plot(label = cluster_types.index, sizes = cluster_types["BGC"], value = cluster_types["BGC"], color=sns.color_palette("muted"), ax = ax, alpha = 0.8, fontsize = 26)
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

def figure_prediction_reactions_violin():
    fn = "../Data/knockouts/all_optknock_results_0.5.csv"
    df = pd.read_csv(fn)
    df["Production increase [%]"] = (df["Production"]-1)*100
    df["Growth"] = -(df["Growth"]-1)*100
    df["Reaction knocked out"] = df["ID"]

    # Sort order
    df_count = df.loc[:, ["ID", "BGC"]].groupby("ID").count()
    df_count.sort_values(by = "BGC", inplace = True, ascending = False)

    fig, ax = plt.subplots(1, figsize = (16,8))
    chart = sns.violinplot(y = "Production increase [%]",  x = "Reaction knocked out", data = df, 
                    ax = ax, scale="count",   cut = 0,  inner="quartile", order = df_count.index)

    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, fontsize = 18)
    sns.despine()
    plt.savefig("../Figures/prediction_reactions_violin.svg")
    

def figure_prediction_reactions_boxplot():
    fn = "../Data/knockouts/all_optknock_results_0.5.csv"
    df = pd.read_csv(fn)
    df["Production increase [%]"] = (df["Production"]-1)*100
    df["Growth"] = -(df["Growth"]-1)*100
    df["Reaction inactivated (by gene knock out)"] = df["ID"]

    # Sort order
    df_count = df.loc[:, ["ID", "BGC"]].groupby("ID").count()
    df_count.sort_values(by = "BGC", inplace = True, ascending = False)

    fig, ax = plt.subplots(1, figsize = (16,8))
    chart = sns.boxplot(y = "Production increase [%]",  x = "Reaction inactivated (by gene knock out)", data = df, 
                        ax = ax, order = df_count.index)

    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, fontsize = 18)
    sns.despine()
    plt.subplots_adjust(bottom = 0.3)
    plt.savefig("../Figures/prediction_reactions_box.svg")
    print(df_count["BGC"])


def prediction_reactions_number():
    fn = "../Data/knockouts/all_optknock_results_0.5.csv"
    df = pd.read_csv(fn)
    df_count = df.loc[:, ["ID", "BGC"]].groupby("ID").count()
    df_count["Number of BGCs"] = df_count["BGC"]
    df_count.reset_index(inplace = True)
    df_count.sort_values(by = "BGC", inplace = True, ascending = False)
    df_count["Reaction inactivated (by gene knock out)"] = df_count["ID"]

    fig, ax = plt.subplots(1, figsize = (16,4))
    sns.barplot(x = "Reaction inactivated (by gene knock out)", y = "Number of BGCs", data = df_count, ax = ax)
    sns.despine()
    plt.subplots_adjust(bottom = 0.3)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, fontsize = 18)
    plt.savefig("../Figures/prediction_reactions_number.svg")





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

    print("#####")
    print(len(taxonomy_list))
    df_tax = pd.DataFrame(taxonomy_list)
    df_tax = df_tax.sort_values(0)
    g1 = df_tax.groupby(1).size()
    #g01 = df_tax.groupby([0, 1]).size().unstack(fill_value = 0)
    g1["Others"] = g1[g1<6].sum()
    g1.name = ""    
    fig, ax = plt.subplots(1, figsize = (10, 10))
    #pie,_,txt = ax.pie(g1[g1>5], radius = 1,autopct='%1.2f%%',  pctdistance=.8)
    squarify.plot(label = g1[g1>5].index, sizes = g1[g1>5], value = g1[g1>5], norm_x = 30, 
                  color = sns.color_palette("muted"), ax = ax, bar_kwargs={'alpha':.8}, text_kwargs={'fontsize':40})
    plt.axis("off")
    plt.savefig("../Figures/taxonomy_treemap.svg")
    g1.to_csv("../Figures/taxonomy_treemap_g.csv")
    df_tax.to_csv("../Figures/taxonomy_treemap_tax.csv")
    df_act = df_tax.loc[df_tax[1]=="Actinobacteria", :]

    g2 = df_act.groupby(2).size()
    g2["Others"] = g2[g2<5].sum()
    g2.name=""
    fig, ax = plt.subplots(1, figsize = (10, 10))
    #pie,_,txt = ax.pie(g1[g1>5], radius = 1,autopct='%1.2f%%',  pctdistance=.8)
    squarify.plot(label = g2[g2>5].index, sizes = g2[g2>5], value = g2[g2>5], norm_x = 30, 
                  color = sns.color_palette("muted"), ax = ax, bar_kwargs={'alpha':.8}, text_kwargs={'fontsize':30})
    plt.axis("off")
    plt.savefig("../Figures/taxonomy_treemap_act.svg")

    df_strep = df_tax.loc[df_tax[3]=="Streptomycetaceae", :]

    g3 = df_strep.groupby(5).size()
    g3["Others"] = g3[g3<5].sum()
    g3.name=""
    fig, ax = plt.subplots(1, figsize = (10, 10))
    #pie,_,txt = ax.pie(g1[g1>5], radius = 1,autopct='%1.2f%%',  pctdistance=.8)
    squarify.plot(label = g3[g3>5].index, sizes = g3[g3>5], value = g3[g3>5], norm_x = 30, 
                  color = sns.color_palette("muted"), ax = ax, bar_kwargs={'alpha':.8}, text_kwargs={'fontsize':30})
    plt.axis("off")
    plt.savefig("../Figures/taxonomy_treemap_strep.svg")

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


def print_names(lst):
    model_fn = "../Models/Sco-GEM.xml"
    model = cobra.io.read_sbml_model(model_fn)
    for r_id in lst:
        r = model.reactions.get_by_id(r_id)
        print(r_id, ";", r.name)


def predict_glycosylation_correlation():
    BGC_count = [5, 3, 4, 1, 2, 2, 1, 1, 2, 2, 2, 5, 2, 3, 2, 2, 1, 1, 1, 1, 1, 2, 1, 5, 2, 1, 1, 2, 2, 1, 1, 1, 3, 0, 5, 3, 2, 2, 6, 3]
    pathway_count = [8, 2, 4, 1, 1, 2, 1, 1, 2, 3, 2, 5, 3, 3, 1, 1, 1, 1, 1, 1, 2, 0, 1, 4, 2, 1, 1, 2, 3, 3, 3, 2, 1, 1, 4, 4, 4, 4, 7, 4]
    (corr, p_value) = pearsonr(pathway_count, BGC_count)
    

    print("Correlation: ", corr)
    print("P-value: ", p_value)
    

if __name__ == '__main__':
    if 0:
        figure_prediction_accuracy_bar_chart()
        figure_prediction_accuracy_full_horizontal()
        figure_prediction_accuracy_full_vertical()
    if 0:
        figure_bigmec_coverage_venn()
        # figure_bigmec_unsuccessful_coverage_venn()
        # figure_bigmec_coverage_pie()
        # figure_bigmec_coverage_treemap()
    if 0:
        figure_prediction_coverage()

    if 0:
        figure_prediction_reactions()
        figure_prediction_reactions_max()
        figure_prediction_reactions_violin()
        figure_prediction_reactions_boxplot()
        prediction_reactions_number()
        lst = ["ASPT", "ASPTA", "GLYCL", "FUM", "MCOALY", "AGT", "GHMT2r", 
              "PGCD", "PSERT", "PSP_L", "TKT1", "ERTHMMOR", "CITMS", "CITCIb", "CITCIa2", "ENO", "PGM"]
        print_names(lst)



    if 0:
        # plot_taxonomic_pie()
        plot_taxonomic_treemap()

    if 0:
        pathway_length_figure()

    if 0:
        figure_production_in_real_vs_constructed_scatter()
    if 0:
        figure_knockouts_predictions_in_real_vs_constructed()
    if 1:
        predict_glycosylation_correlation()
