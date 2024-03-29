#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fileencoding=utf-8

"""
TRANSLATE FASTA
James B. Pease
"""

import os
import sys
import argparse
from itertools import groupby
from random import sample as rsample
from Bio.SubsMat import MatrixInfo
# import numpy as np
import time

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

AMBIG = {"A": "A",
         "C": "C",
         "G": "G",
         "T": "T",
         "R": "AG",
         "Y": "CT",
         "W": "AT",
         "S": "CG",
         "M": "AC",
         "K": "GT",
         "B": "CGT",
         "D": "AGT",
         "H": "ACT",
         "V": "ACG",
         "N": "ACGT"}


def fasta_iter(fasta_name):
    """
        given a fasta file. yield tuples of header, sequence
        Adapted from https://github.com/brentp
    """
    filehandler = open(fasta_name, 'r')
    faiter = (x[1] for x in groupby(filehandler, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


def calc_pwdist_nuc(seq0, seq1, mode="random2", alpha=1,
                    minoverlap=200):
    diff = 0
    total = 0
    for i, base0 in enumerate(seq0):
        base0 = base0.upper()
        base1 = seq1[i].upper()
        if (base0 not in 'RYWSMKBDHVNATGC' or
                base1 not in 'RYWSMKBDHVNATGC'):
            continue
        if base0 in "BDHVN" or base1 in "BDHVN":
            if mode != "random4":
                continue
        if base0 in "RYWSMK" or base1 in "RYWSMK":
            if mode != "random2":
                continue
        if mode in ("random2", "random4"):
            base0 = rsample(AMBIG[base0], 1)[0]
            base1 = rsample(AMBIG[base1], 1)[0]
        total += 1
        if base0 + base1 in (
            'AC', 'AT', 'CA', 'CG',
                'GC', 'GT', 'TA', 'TG'):
            diff += alpha
        elif base0 + base1 in (
                'AG', 'CT', 'GA', 'TC'):
            diff += 1
    if total < minoverlap:
        return "na"
    return float(diff / total)


def calc_pwdist_prot(seq0, seq1, model="simple",
                     minoverlap=200):
    diff = 0
    total = 0
    matrix = None
    if model.startswith('blosum'):
        if model == "blosum62":
            matrix = MatrixInfo.blosum62
        elif model == "blosum45":
            matrix = MatrixInfo.blosum45
        elif model == "blosum80":
            matrix = MatrixInfo.blosum80
        for i, base0 in enumerate(seq0):
            base0 = base0.upper()
            base1 = seq1[i].upper()
            if (base0 not in 'ACDEFGHIKLMNPQRSTVWY' or
                    base1 not in 'ACDEFGHIKLMNPQRSTVWY'):
                continue
            if base0 != base1:
                pair = (base0, base1)
                if pair in matrix:
                    diff += matrix[pair]
                else:
                    diff += matrix[(tuple(reversed(pair)))]
                total += 1
    else:
        for i, base0 in enumerate(seq0):
            base0 = base0.upper()
            base1 = seq1[i].upper()
            if (base0 not in 'ACDEFGHIKLMNPQRSTVWY' or
                    base1 not in 'ACDEFGHIKLMNPQRSTVWY'):
                continue
            total += 1
            if base0 != base1:
                diff += 1
    if total < minoverlap:
        return "na"
    return float(diff / total)


def generate_argparser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=_LICENSE)
    parser.add_argument('fasta', type=os.path.abspath,
                        help="input fasta file")
    parser.add_argument('--out', type=os.path.abspath, required=True,
                        help="output fasta file")
    parser.add_argument("--mode",
                        choices=("strict", "random2", "random4"),
                        default="random2",
                        help=("strict=skip ambiguous; "
                              "random2=choose random at biallelic only; "
                              "random4=choose random at any ambiguous"))
    parser.add_argument("--seqtype", choices=("nuc", "prot"),
                        default="nuc",
                        help="nucleotide or protein alignment")
    parser.add_argument("--min-overlap", type=int, default=200)
    parser.add_argument("--model", choices=("simple",
                                            "blosum62", "blosum45",
                                            "blosum80"),
                        default="simple",
                        help="simple distance or amino acid model")
    parser.add_argument("--alpha", default=1,
                        help=("transversion/transition score ratio, change"
                              " from alpha=1 to set K80 model"))
    return parser


def main(arguments=None):
    time0 = time.time()
    arguments = arguments if arguments is not None else sys.argv[1:]
    parser = generate_argparser()
    args = parser.parse_args(args=arguments)
    headers = []
    seqs = {}
    for hdr, seq in fasta_iter(args.fasta):
        headers.append(hdr)
        seqs[hdr] = seq
    with open(args.out, 'w') as outfile:
        for i, hdr0 in enumerate(headers):
            for j in range(i+1, len(headers)):
                hdr1 = headers[j]
                if args.seqtype == "nuc":
                    dist = calc_pwdist_nuc(
                        seqs[hdr0], seqs[hdr1],
                        mode=args.mode,
                        alpha=args.alpha,
                        minoverlap=args.min_overlap)
                elif args.seqtype == "prot":
                    dist = calc_pwdist_prot(
                        seqs[hdr0], seqs[hdr1],
                        model=args.model,
                        minoverlap=args.min_overlap)

                outfile.write("\t".join([hdr0, hdr1, str(dist)]) + "\n")
    print(time.time() - time0)
    return ''


if __name__ == '__main__':
    main()
