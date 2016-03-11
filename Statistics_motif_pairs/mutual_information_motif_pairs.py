import csv
from collections import OrderedDict
import numpy as np
import os

def arguments():
    """adding arguments for the script"""
    import argparse
    parser = argparse.ArgumentParser(description="""
    Given a motif name, it computes the mutual information between the
    motif and the other motifs. Optionally, it receives a set of regions
    to limit the calculation of mutual information only for those subsets.
    """)
    parser.add_argument('-o', '--ovr', dest='overlapping_dir', action='store',
                        type=str, required=True,
                        help="""A directory with results for counting
                        overlapping motifs. Made by ''count_overlapping_motifs.py""")
    parser.add_argument('-n', '--nvr', dest='nonoverlapping_dir', action='store',
                        type=str, required=True,
                        help="""A directory with results for counting
                        non-overlapping motifs. Made by ''count_nonoverlapping_motifs.py""")
    parser.add_argument('-i', '--motevo', dest='motevo_dir', action='store',
                        type=str, required=True,
                        help="""Directory with all MotEvo output files""")
    parser.add_argument('-m', '--motif', dest='motif', action='store',
                        type=str, required=True,
                        help="""Motif name, it's bascially name of one of the directory in
                        the overlapping or non-overlapping directory.""")
    parser.add_argument('-p', '--proxy', dest='proxyBED', action='store',
                        type=str, required=False,
                        help="""Optional BED file that contains regions of interest.
                        Mutual information is only calculated for the these regions.""")
    args = parser.parse_args()
    return args


def main():
    
    args = arguments()
    if args.proxy:
        regions_of_interest = dict([(l.split()[4], 0)
                                    for l in open(args.proxy)])   # by default it' assumed that 5th
    else:
        regions_of_interest = None

    motif_names = os.listdir(args.motevo_dir)
    print motif_names


if __name__ == '__main__':
    main()


