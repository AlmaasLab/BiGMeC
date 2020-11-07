#!/usr/bin/env python
# coding: utf-8
"""
Authors:
  - Fredrik Fossheim
  - Snorre Sulheim, snorres.sulheim@sintef.no


Date: 24.09.2020
Lisence: CC-BY-4.0

This file is used to get the gbk-files from antismash
"""

from Bio import SeqIO
import re
import os
import json
import test
import csv
from bs4 import BeautifulSoup
from urllib.request import urlopen
from pathlib import Path


def get_web_data(url):
    html = urlopen(url)
    soup = BeautifulSoup(html, features="lxml")
    return soup.body.p

def write_gbk_file(gbk_entry, path):
    with open(path, "w+") as f:
        f.write(gbk_entry)
    f.close()

def get_all(folder):

    folder = Path(folder)
    folder.mkdir(exist_ok = True)

    for a in range(1,2070,1):
        # path er bare pathen til en tom mappe som alle GBK-filene lagres i
        path = str(folder / "BGC{0}.gbk".format(str(a).zfill(7)))
        print(path)
        realUrl = 'https://mibig.secondarymetabolites.org/repository/BGC' + str(a).zfill(7) + '/generated/BGC' \
                  + str(a).zfill(7) + '.1.region001.gbk'

        # putter selve datainnhentinga i en try, fordi det er en del av linkene som ikke finnes, og da crasher programmet
        # når man kommer til GBK nummer 15, som ikke finnes i MiBiG
        try:
           string = str(get_web_data(realUrl))[3:-4]#this is the entire gbk-file for each cluster
        except Exception:
            print(a)
            continue
        # Write GBK file to specified folder path
        write_gbk_file(string, path)

if __name__ == '__main__':
    folder = "../Data/MiBiG2"
    get_all(folder)