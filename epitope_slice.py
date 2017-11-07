#!/usr/bin/env python

"""
Parses output from vep_parse.py to produce epitopes.
"""

# standard libs
import argparse
import gzip
import re
import sys
# third party libs
from Bio import SeqIO
import Bio.Data.CodonTable


__author__ = 'dkdeconti'
__copyright__ = "Copyright 2017"
__credits__ = ["Derrick DeConti"]
__license__ = "MIT"
__maintainer__ = "Derrick DeConti"
__email__ = "deconti@jimmy.harvard.edu"
__status__ = "Development"


def collapse_epitopes(epitopes):
    """
    Inverts the epitope map and merges epitopes with enst and pos added.
    """
    output_epitopes = {}
    for k, v in epitopes.items():
        enst, ensg, name, pos = k
        #alt_epitope, wt_epitope = v
        if v in output_epitopes:
            output_epitopes[v][1].append((enst, str(pos)))
        else:
            output_epitopes[v] = [name, [(enst, str(pos))], ensg]
    return output_epitopes


def get_epitope(id_pos_map, protein_map, flank):
    """
    Retrieves a FASTA of protein seq from list.
    """
    epitopes = {}
    for k, v in id_pos_map.items():
        try:
            protein_seq = protein_map[k]
        except KeyError:
            continue
        enst = k
        ensg, name, pos, aa_change = v
        begin = int(pos) - flank
        if begin < 0:
            begin = 0
        end = int(pos) + flank + 1
        alt_epitope = ''.join([protein_seq[begin:pos],
                               aa_change,
                               protein_seq[pos+1:end]])
        wt_epitope = protein_seq[begin: end]
        epitopes[(enst, ensg, name, pos)] = (alt_epitope, wt_epitope)
    return epitopes


def parse_input_tsv(filename):
    """
    Parses the output from vep_parse to output map of id and pos.
    """
    output_map = {}
    with open(filename, 'rU') as handle:
        for line in handle:
            arow = line.strip('\n').split('\t')
            if len(arow) <= 1: continue # skip empty lines
            # skips if insertion or stop mutation
            if len(arow[4]) == 3 or arow[4][-1] != '*':
                output_map[arow[1]] = (arow[0], arow[2], int(arow[3]), arow[4][-1])
    return output_map


def parse_cds_fasta(filename):
    """
    Parses the Ensembl CDS FASTA and translates CDS to amino acids.
    """
    output_map = {}
    # No "with" because biopython chokes on it
    handle = open(filename, 'rU')
    for seq_record in SeqIO.parse(handle, "fasta"):
        key = seq_record.id.split('.')[0]
        try:
            value = str(seq_record.seq.translate(cds=True))
            output_map[key] = value
        except Bio.Data.CodonTable.TranslationError as err:
            #print key, err
            pass
    handle.close()
    return output_map


def write_epitopes(epitopes, alt_filename, wt_filename):
    """
    Writes epitopes as FASTA to alt and wt files.
    """
    with open(alt_filename, 'w') as out_alt, open(wt_filename, 'w') as out_wt:
        for k, v in epitopes.items():
            alt_epitope, wt_epitope = k
            gene, enst_list, ensg = v
            gene_name = '-'.join([gene,])
            enst = [' : '.join(i) for i in enst_list]
            desc = ' | '.join([gene,
                               ' , '.join(enst),
                               ensg])
            out_alt.write('>' + desc + '\n' + alt_epitope + '\n')
            out_wt.write('>' + desc + '\n' + wt_epitope + '\n')


def main():
    """
    Arg parsing and central dispatch.
    """
    # Arg parsing
    desc = "Parses output from vep_parse.py to produce epitopes"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--cds", metavar="CDS_FASTA", required=True,
                        help="NCBI CDS FASTA file")
    parser.add_argument("--alt", metavar="OUTPUT_ALT", required=True,
                        help="filename for alt epitope output FASTA")
    parser.add_argument("--wt", metavar="OUTPUT_wt", required=True,
                        help="filename for wt epitope output FASTA")
    parser.add_argument("-f", "--flanking", metavar="BP",
                        default=8, type=int,
                        help="flanking +/- bp")
    parser.add_argument("input", metavar="INPUT")
    args = parser.parse_args()
    # Central dispatch
    id_pos_map = parse_input_tsv(args.input)
    protein_map = parse_cds_fasta(args.cds)
    epitopes = collapse_epitopes(get_epitope(id_pos_map, protein_map,
                                             args.flanking))
    write_epitopes(epitopes, args.alt, args.wt)
    


if __name__ == "__main__":
    main()
