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
# plt.style.use("seaborn")

plt.rcParams.update({'font.size': 22})

prediction_scores = {
    "Domain function \nin total": [177, 215],
    "Extender units": [72, 92],
    "Non-extending \ndomains": [103, 123],
    "Starter units": [5, 8]}

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

def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

if __name__ == '__main__':
    figure_prediction_accuracy_bar_chart()
