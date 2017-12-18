#! /usr/bin/env python

"""
Parses the NetHMC output into TSV sorted by affinity.
"""

import argparse
import itertools
import re
import sys
from operator import itemgetter


def parse_fasta_for_names(filename):
    """

    """
    name_dict = {}
    with open(filename, 'rU') as handle:
        for line in handle:
            if not re.match(">", line):
                continue
            name = line.strip('\n')[1:]
            name_sub = name[:15]
            name_dict[name_sub] = name
    return name_dict


def parse_tsv(filename, name_dict):
    """

    """
    output_matrix = []
    with open(filename, 'rU') as handle:
        curr_protein = []
        for line in handle:
            if line[0] == "#" or line[0] == "-" or len(line.strip('\n')) < 1:
                continue
            if re.match("Protein", line):
                continue
            arow = line.strip('\n').split()
            if arow[0] == "pos":
                continue
            arow[12] = float(arow[12])
            if len(arow[10].split('-')) == 3:
                #arow = arow[:10] + arow[10].split('_') + arow[11:]
                arow = arow[:10] + name_dict[arow[10]].split('-') + arow[11:]
                #print arow
            output_matrix.append(arow)
    return output_matrix


def main():
    """
    Arg parsing and central dispatch.
    """
    # arg parsing
    parser = argparse.ArgumentParser(description="Parse NetHMC output to TSV.")
    parser.add_argument("-f", "--fasta", metavar="FASTA",
                        help="fasta file for epitopes input into nethmc")
    parser.add_argument("nethmc", metavar="NETHMC_OUTPUT",
                        nargs='+',
                        help="Output file from NetHMC")
    args = parser.parse_args()
    # Central dispatch
    name_dict = parse_fasta_for_names(args.fasta)
    lom = list(itertools.chain.from_iterable([parse_tsv(nethmc, name_dict)
                                              for nethmc in args.nethmc]))
    #print set([len(v) for v in list(itertools.chain.from_iterable(lom))])
    nethmc_matrix = sorted([v for v in lom if len(v) >= 18 and v[16] <= 50],
                           key=itemgetter(14))
    sys.stdout.write('\t'.join(["pos", "HLA", "peptide", "Core Offset", "I_pos",
                                "I_len", "D_pos", "D_len", "iCore",
                                "identity", "chrom", "genomic_pos",
                                "gene_symbol", "ref_allele", "alt_allele", 
                                "1-log50k(aff)", "Affinity(nM)",
                                "%Rank BindLevel"]) + '\n')
    for v in nethmc_matrix:
        sys.stdout.write('\t'.join([str(i) for i in v[:18]]) + '\n')
    #nethmc_matrix


if __name__ == "__main__":
    main()
