#!/usr/bin/env python

"""
Parses output from vep_parse.py to produce epitopes.
"""

import argparse
import gzip
import re
import sys


__author__ = 'dkdeconti'
__copyright__ = "Copyright 2017"
__credits__ = ["Derrick DeConti"]
__license__ = "MIT"
__maintainer__ = "Derrick DeConti"
__email__ = "deconti@jimmy.harvard.edu"
__status__ = "Development"

''' Example code to reference
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "me@mysite.com"
id_list = set(open('test.txt', 'rU'))
handle = Entrez.efetch(db="protein", id=id_list, rettype="fasta", retmode="text")   
for seq_record in SeqIO.parse(handle, "fasta"):
    print ">" + seq_record.id, seq_record.description
print seq_record.seq
handle.close()
'''

def main():
    """
    Arg parsing and central dispatch.
    """
    # Arg parsing
    desc = "Parses output from vep_parse.py to produce epitopes"
    parser = argparse.ArgumentParser(description=desc)

    argse = parser.parse_args()
    # Central dispatch
    pass


if __name__ == "__main__":
    main()
