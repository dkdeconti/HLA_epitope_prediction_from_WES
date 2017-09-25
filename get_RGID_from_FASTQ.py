#!/usr/bin/python

"""
Pull RGID information (from CCCB NextSeq naming scheme) from FASTQ.
"""

import argparse
import gzip


def get_rg_values(filename, sample_name):
    """
    Parses gzipped fastq for RG values and returns RG str for BWA.
    """
    with gzip.open(filename, 'rU') as handle:
        line = handle.readline()
        arow = line.strip('\n').split()
        info = arow[0].split(':')[1:]
        instrument_id = info[0]
        run_id = info[1]
        flowcell_id = info[2]
        flowcell_lane = info[3]
        index_seq = arow[1].split(':')[3]
        rgid = '.'.join([sample_name, flowcell_id, flowcell_lane])
    rglb = '.'.join([sample_name, run_id])
    rgpu = '.'.join([instrument_id,
                     flowcell_lane,
                     index_seq])
    rgsm = sample_name
    rgcn = "DFCI-CCCB"
    rgpl = "ILLUMINA"
    rg_vals = "@RG\\tID:" + rgid + "\\tPL:" + rgpl + "\\tLB:" + \
              rglb + "\\tSM:" + rgsm + "\\tCN:" + rgcn + "\\tPU:" + rgpu
    return rg_vals


def main():
    """
    Arg parsing and central dispatch.
    """
    # Arg parsing
    desc = "Pull RGID info from CCCB NextSeq produced FASTQ"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("samplename", metavar="SAMPLENAME", help="sample name")
    parser.add_argument("fastq", metavar="FASTQ", help="input fastq file")
    args = parser.parse_args()
    # Central dispatch
    print get_rg_values(args.fastq, args.samplename)


if __name__ == "__main__":
    main()
