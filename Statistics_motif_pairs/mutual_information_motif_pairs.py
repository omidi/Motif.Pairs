import csv
from collections import OrderedDict
import numpy as np
from numpy import log, exp
import os
from scipy.special import gammaln
import re


def arguments():
    """adding arguments for the script"""
    import argparse
    parser = argparse.ArgumentParser(description="""
    Given a motif name, it computes the mutual information between the
    motif and the other motifs. Optionally, it receives a set of regions
    to limit the calculation of mutual information only for those subsets.
    """)
    parser.add_argument('-m', '--motif', dest='motif', action='store',
                        type=str, required=True,
                        help="""Motif name, it's bascially name of one of the directory in
                        the overlapping or non-overlapping directory.""")
    # parser.add_argument('-n', '--novr', dest='nonoverlapping_dir', action='store',
    #                     type=str, required=True,
    #                     help="""A directory with results for counting
    #                     non-overlapping motifs. Made by ''count_nonoverlapping_motifs.py""")
    # parser.add_argument('-i', '--motevo', dest='motevo_dir', action='store',
    #                     type=str, required=True,
    #                     help="""Directory with all MotEvo output files""")
    parser.add_argument('-r', '--regions', dest='regions', action='store',
                        type=str, required=True,
                        help="""A BED file that contains regions of interest.
                        Mutual information is only calculated for the these regions.""")
    parser.add_argument('-c', '--cutoff', dest='cutoff', action='store',
                        type=float, required=False, default=0.5,
                        help="""Posterior cutoff for selecting TFBS from the MotEvo output.
                        (Default 0.5)""")
    args = parser.parse_args()
    if args.cutoff>1 or args.cutoff<0:
        print 'Posterior cutoff is invald!'
        print 'Program halts!'
        exit()
    return args


def load_double_counts_matrix(dirname, motif, regions):
    fname = os.path.join(dirname, motif, motif)
    with open(fname) as inf:
        header = dict([(m, i) for i, m in enumerate(inf.readline().split())])
        double_counts = dict([(r, np.zeros(len(header), dtype=np.int)) for r in regions])
        for line in inf:
            region = line.split()[0]
            if region in double_counts:
                double_counts[region] = map(int, line.split()[1:])
    double_counts = np.transpose(np.array([double_counts[r] for r in regions]))  # order according to the orders in regions and make matrix
    return header, double_counts


def load_single_counts(dirname, motif, regions, cutoff=.5):
    """
    receieves the MotEvo directory and a name of a motif and
    a set of regions. Then by setting a posterior cutoff
    of 0.5 it extracts the vector of single counts.
    """
    fname = os.path.join(dirname, motif)
    counts = dict([(r, 0) for r in regions])
    with open(fname) as inf:
        for line in inf:
            rec = line.split()
            region = rec[3].split(";")[-1].rstrip()
            if region in counts:
                if float(rec[4]) > cutoff:
                    counts[region] = 1
    return np.array([counts[r] for r in regions])


def correlation_score(N1, N2):
    pseudo = .5
    n1 = np.sum(N1)
    n2 = np.sum(N2)
    score = log(2.) + log(pseudo) + gammaln(len(N1) + 2*pseudo) + gammaln(len(N2) + 2*pseudo)
    score += gammaln(n1 + n2 + pseudo) + gammaln(2*len(N1) - n1 - n2 + pseudo)
    score -= log(2.*pseudo) + gammaln(2*len(N1) + pseudo*2) + gammaln(n1 + pseudo) + gammaln(n2 + pseudo)
    score -= gammaln(len(N1) - n1 + pseudo) + gammaln(len(N2) - n2 + pseudo)
    return score


def mutual_information(N1, N2):
    pseudo = .5
    p1 = {1:(np.sum(N1)+pseudo)/(len(N1) + 2*pseudo),
          0:(len(N1) - np.sum(N1)+pseudo)/(len(N1) + 2*pseudo)}
    p2 = {1:(np.sum(N2)+pseudo)/(len(N2) + 2*pseudo),
          0:(len(N2) - np.sum(N2)+pseudo)/(len(N2) + 2*pseudo)}
    p12 = {(1,1):(len([1 for n1, n2 in zip(N1,N2) if n1==1 and n2==1]) + pseudo) / (len(N1) + 4*pseudo),
           (1,0):(len([1 for n1, n2 in zip(N1,N2) if n1==1 and n2==0]) + pseudo) / (len(N1) + 4*pseudo),
           (0,1):(len([1 for n1, n2 in zip(N1,N2) if n1==0 and n2==1]) + pseudo) / (len(N1) + 4*pseudo),
           (0,0):(len([1 for n1, n2 in zip(N1,N2) if n1==0 and n2==0]) + pseudo) / (len(N1) + 4*pseudo),
           }
    I = 0
    for n1, n2 in zip(N1, N2):
        I += p12[(n1, n2)] * np.log2( p12[(n1, n2)] / (p1[n1]*p2[n2]) )
    return I



def main():
    args = arguments()
    args.motevo_dir = "/home/somidi/scratch/Tissue.Specificity/Motif.Pairs/MotEvo.Outputs/min.post.50/"
    args.overlapping = "/home/somidi/scratch/Tissue.Specificity/Motif.Pairs/Pairs.Counts/Overlapping/"
    regions = [l.split()[4] for l in open(args.regions)]   # by default it's assumed that 5th contains region IDs
    motif_names = os.listdir(args.motevo_dir)
    header, double_counts = \
            load_double_counts_matrix(args.overlapping, args.motif, regions)
    # motif_names = [m for m in motif_names if re.search('RXR', m) ]
    N1 = load_single_counts(args.motevo_dir, args.motif, regions, args.cutoff)
    for motif in motif_names:
        N2 = load_single_counts(args.motevo_dir, motif, regions, args.cutoff)
        print '\t'.join([motif,
                         "%0.4f" % mutual_information(N1, N2),
                         "%d" % np.sum(double_counts[header[motif], ...]),
                         "%d" % len(regions),
                         ])


if __name__ == '__main__':
    main()


