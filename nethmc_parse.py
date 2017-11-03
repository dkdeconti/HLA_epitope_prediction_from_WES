#! /usr/bin/env python

"""
Parses the NetHMC output into TSV sorted by affinity.
"""

import argparse
import itertools
import re
import sys
from operator import itemgetter


def parse_tsv(filename):
    """

    """
    output_matrix = []
    with open(filename, 'rU') as handle:
        for line in handle:
            if line[0] == "#" or line[0] == "-" or len(line.strip('\n')) < 1:
                continue
            if re.match("Protein", line):
                continue
            arow = line.strip('\n').split()
            if arow[0] == "pos":
                continue
            #print line.strip('\n').split()
            arow[12] = float(arow[12])
            output_matrix.append(arow)
    return output_matrix


def main():
    """
    Arg parsing and central dispatch.
    """
    # arg parsing
    parser = argparse.ArgumentParser(description="Parse NetHMC output to TSV.")
    parser.add_argument("fasta", metavar="FASTA", required=True,
                        help="fasta file for epitopes input into nethmc")
    parser.add_argument("nethmc", metavar="NETHMC_OUTPUT",
                        nargs='+',
                        help="Output file from NetHMC")
    args = parser.parse_args()
    # Central dispatch
    lom = [parse_tsv(nethmc) for nethmc in args.nethmc]
    nethmc_matrix = sorted(list(itertools.chain.from_iterable(lom)),
                           key=itemgetter(12),
                           reverse=False)
    sys.stdout.write('\t'.join(["pos", "HLA", "peptide", "Core Offset", "I_pos",
                                "I_len", "D_pos", "D_len", "iCore",
                                "identity 1-log50k(aff)", "Affinity(nM)",
                                "%Rank BindLevel"]) + '\n')
    for v in nethmc_matrix:
        sys.stdout.write('\t'.join([str(i) for i in v]) + '\n')
    #nethmc_matrix



if __name__ == "__main__":
    main()
