#!/usr/bin/env python
"""
Parses VEP annotated VCF for rare variants by ExAC.
"""

import argparse
import gzip
import re
import sys


__author__ = 'konradjk, dkdeconti'
__copyright__ = "Copyright 2017"
__credits__ = ["Konrad Karczewski", "Derrick DeConti"]
__license__ = "MIT"
__maintainer__ = "Derrick DeConti"
__email__ = "deconti@jimmy.harvard.edu"
__status__ = "Development"


def main(args):
    f = gzip.open(args.vcf) if args.vcf.endswith('.gz') else open(args.vcf)
    field_names = None
    header = None
    for line in f:
        line = line.strip()

        # Reading header lines to get VEP and individual arrays
        if line.startswith('#'):
            print line
            line = line.lstrip('#')
            #print line
            if line.find('ID=CSQ') > -1:
                field_names = line.split('Format: ')[-1].strip('">').split('|')
            if line.startswith('CHROM'):
                header = line.split()
                # Creates dict of index and column header
                header = dict(zip(header, range(len(header))))
                #print header
            continue

        if field_names is None:
            err_output = "VCF file does not have a VEP header line. Exiting."
            print >> sys.stderr, err_output
            sys.exit(1)
        if header is None:
            err_output = "VCF file does not have a header line" + \
            "(CHROM POS etc.). Exiting."
            print >> sys.stderr, err_output
            sys.exit(1)

        # Pull out annotation info from INFO and ALT fields
        fields = line.split('\t')
        info_field = dict([(x.split('=', 1))
                            if '=' in x else (x, x)
                            for x in re.split(';(?=\w)',
                                              fields[header['INFO']])])
        print fields
        # Only reading lines with an annotation after this point
        if 'CSQ' not in info_field: continue
        annotations = [dict(zip(field_names, x.split('|')))
                       for x in info_field['CSQ'].split(',')
                       if len(field_names) == len(x.split('|'))]
        consequences = [annot["Consequence"] for annot in annotations]
        if 'missense_variant' in consequences:
            print consequences
        exac_freq = [float(freq) if freq != '' else 0.0
                     for freq in info_field['CSQ'].split('|')[-9:-5]]        
        if max(exac_freq) == 0.0:
            for x in annotations:
                if x['Amino_acids'] != "":
                    print '\t'.join([x['Gene'], x['Feature'], x['SYMBOL'],
                                     x['Protein_position'], x['Amino_acids']])
        # Explanations #######################################################
        # annotations = list of dicts, each corresponding to a 
        # transcript-allele pair
        # (each dict in annotations contains keys from field_names)
        # lof_annotations = list of dicts, only high-confidence LoF
        # if len(lof_annotations) > 0 can determine if at least one allele is 
        # LoF for at least one transcript
        # fields = vcf line, can be accessed using:
        # fields[header['CHROM']] for chromosome,
        # fields[header['ALT']] for alt allele,
        # or samples using sample names, as fields[header['sample1_name']]
        ######################################################################
    f.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf', '--input', '-i',
                        help='Input VCF file (from VEP+LoF); may be gzipped',
                        required=True)
    parser.add_argument('-m', '--minimum', metavar="PVALUE",
                        type=float, default=0.05,
                        help="Minimum ExAC_all freq for pass")
    args = parser.parse_args()
    main(args)