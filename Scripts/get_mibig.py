#!/usr/bin/env python
# coding: utf-8
"""
Copyright 2020 Snorre Sulheim (snorre.sulheim@sintef.no)
https://github.com/AlmaasLab/BiGMeC

This file is used to extract the antismash results for all BGCs at MiBiG (up to BGC0002070).


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
    """
    Get the antiSMASH results for all BGCs (up to BGC0002070) from MIBiG.
    Store the results in the given folder.
    """
    folder = Path(folder)
    folder.mkdir(exist_ok = True)

    for a in range(1,2070,1):
        path = str(folder / "BGC{0}.gbk".format(str(a).zfill(7)))
        print(path)
        realUrl = 'https://mibig.secondarymetabolites.org/repository/BGC' + str(a).zfill(7) + '/generated/BGC' \
                  + str(a).zfill(7) + '.1.region001.gbk'

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