#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fileencoding=utf-8

"""
Plot an MDS from Pairwise Matrix
James B. Pease
"""

import os
import sys
import argparse
import numpy as np
import matplotlib.colors
from matplotlib.colors import ListedColormap
from sklearn.manifold import MDS
from matplotlib import pyplot as plt, cm, colors
from time import time

_LICENSE = """
http://www.github.org/jbpease/mixtape
MixTAPE: Mix of Tools for Analysis in Phylogenetics and Evolution

This file is part of MixTAPE.

MixTAPE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MixTAPE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MixTAPE.  If not, see <http://www.gnu.org/licenses/>.
"""

np.seterr(divide='ignore', invalid='ignore')

class CSVfile(object):
    "field-delimited data file"

    def __init__(self, fpath, delimiter=",",
                 colheaders=True):
        self.path = os.path.abspath(fpath)
        self.delimiter = delimiter
        self.colheaders = colheaders
        if self.colheaders is True:
            with open(self.path) as csvfile:
                self.colheaders = csvfile.readline().split(
                    self.delimiter)

    def __iter__(self):
        firstline = True
        with open(self.path) as csvfile:
            for line in csvfile:
                if self.colheaders is not False and firstline is True:
                    firstline = False
                    continue
                yield line.rstrip().split(self.delimiter)



COLORS = {'P1': 'xkcd:violet',
          'P2': 'g',
          'P3': 'b',
          'P4': 'orange',
          'P6': 'm',
          'P17': 'r',
          'P30': 'r',
          'P31': 'r',
          'P35': 'r',
          'P36': 'r',
          'P40': 'r',
          'P42': 'xkcd:bubblegum',
          'P43': 'xkcd:bubblegum',
          'P43': 'xkcd:bubblegum',

         }
def generate_argparser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=_LICENSE)
    parser.add_argument('pwdist', type=os.path.abspath,
                        help="input pairwise list")
    parser.add_argument("--csv", action="store_true",
                        help="comma-delimited (default:tab-delimited)")
    parser.add_argument("--labeldelim", default="/",
                        help="label delimeter")
    parser.add_argument("--labelfield", default=6, type=int,
                        help="label delimited fields")
#    parser.add_argument("--groupfield", type=int,
#                        help="label delimited field for the group for the key")
    parser.add_argument('--out', type=os.path.abspath,
                        default=sys.stdout, help="output fasta file")
    parser.add_argument('--figsize', type=float, nargs=2,
                        help="fixed height/width in inches for figure")
    return parser


def main(arguments=None):
    arguments = arguments if arguments is not None else sys.argv[1:]
    parser = generate_argparser()
    args = parser.parse_args(args=arguments)
    time0 = time()
    #nseqs = len(csvfile.colheaders) - 1
    #arr = np.zeros((nseqs, nseqs))
    #cmap = cm.tab20
    colors = ('xkcd:puke green', 'xkcd:sky blue', 'xkcd:gold', 'xkcd:rust', 'xkcd:light violet', 'xkcd:chartreuse', 'xkcd:cyan', 'xkcd:light tan', 'xkcd:coral', 'xkcd:bluish purple', 'xkcd:grass', 'xkcd:turquoise', 'xkcd:bright yellow', 'xkcd:maroon', 'xkcd:grape', 'xkcd:jungle green', 'xkcd:navy blue', 'xkcd:golden yellow', 'xkcd:dark pink', 'xkcd:periwinkle blue', 'xkcd:dark olive green', 'xkcd:royal blue', 'xkcd:pastel orange', 'xkcd:light pink', 'xkcd:medium purple', 'xkcd:poop', 'xkcd:dusty rose', 'xkcd:hospital green', 'xkcd:bright red', 'xkcd:faded yellow', 'xkcd:light aqua', 'xkcd:diarrhea', 'xkcd:cerise', 'xkcd:melon', 'xkcd:bubblegum')
    cmap = matplotlib.colors.ListedColormap(colors)
    norm = matplotlib.colors.Normalize(vmin=1, vmax=35, clip=False)
    #class matplotlib.colors.ListedColormap(colors,name='from_List',N=35)
    #cmap = cm.from_List
    if args.figsize is not None:
        fig = plt.figure(figsize=(args.figsize[0], args.figsize[1]), )
    else:
        fig = plt.figure()
    ax1 = fig.add_subplot(111)
    i = 0
    headers = []
    first_entry = None
    nseq = 0
    headers = []
    hindex = {}
    with open(args.pwdist) as infile:
        for line in infile:
            entry = line.rstrip().split(
                "\t" if args.csv is False else ",")
            if first_entry is None:
                first_entry = entry[0]
                headers.append(entry[0])
                hindex[entry[0]] = 0
                headers.append(entry[1])
                hindex[entry[1]] = 1
                nseq += 2
                continue
            if entry[0] != first_entry:
                break
            headers.append(entry[1])
            hindex[entry[1]] = nseq - 1
            nseq += 1
    print(nseq)
    seed = np.random.RandomState(seed=3)
    arr = np.zeros((nseq, nseq))
    with open(args.pwdist) as infile:
        for line in infile:
            entry = line.rstrip().split()
            if entry[2] == "na":
                continue
            i = hindex[entry[0]]
            j = hindex[entry[1]]
            arr[i][j] = float(entry[2])
            arr[j][i] = float(entry[2])
    #for i in range(nseq):
    #    arr[i][i] = 1.0

    mds = MDS(n_components=2, dissimilarity='precomputed',
              max_iter=3000,
              random_state=seed,
              eps=1e-9)
    #arr = np.exp(arr)
    print(arr)
    results = mds.fit(arr)
    coords = results.embedding_
    #if args.groupfield is not None:
        #global rep
    #    rep = set([])
    #for label, x, y in zip(csvfile.colheaders[1:], coords[:, 0], coords[:, 1]):
    grouplabels = [x.split(args.labeldelim)[args.labelfield] for x in headers]
    color_index = {}
    icolor = 0
    for x in grouplabels:
        if x not in color_index:
            color_index[x] = colors[icolor]
            icolor += 1
            if icolor > 34:
                icolor = 34

    #groupcolors = [color_index[x] for x in grouplabels]
    groupcolors = [COLORS.get(x, 'k') for x in grouplabels]
    ax1.scatter(coords[:, 0], coords[:, 1], marker='o', s=3, c=groupcolors)

    for i, x in enumerate(grouplabels):
        ax1.annotate(x, xy=(coords[i, 0], coords[i, 1]), size="4")

     #       print(x, y)
     #       ax1.annotate(
     #           #grplabel,
     #           label[1], xy=(x, y), xytext=(-1, 1), size="4",
     #               textcoords='offset points', ha='right', va='center')
     #   if grplabel not in rep:
     #       rep.update([grplabel])
    #plt.show()
    #ax1.legend(prop={'size': 10})
    plt.tight_layout()
    plt.savefig(args.out)
    print('Completed in {:.2f} seconds.'.format(time() - time0))
    return ''


if __name__ == '__main__':
    main()

