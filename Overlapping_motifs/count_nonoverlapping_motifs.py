import csv
from collections import OrderedDict
import numpy as np
import os

def arguments():
    """adding arguments for the script"""
    import argparse
    parser = argparse.ArgumentParser(description="""
    Given a directory of MotEvo output files, and a filename
    that represents a motif, it calculates all overlapping
    motif counts for the given motif versus all the other motifs
    in the directory (including the same motif with itself).
    """)
    parser.add_argument('-d', '--dir', dest='dirname', action='store',
                        type=str, required=True,
                        help="""MotEvo output directory, containing all the output
                        files from MotEvo. """)
    parser.add_argument('-i', '--input', dest='input_file', action='store',
                        type=str, required=True,
                        help='The filename for the motif of interest in the '
                             'MotEvo output directory')
    parser.add_argument('-o', '--output', dest='output_dir', action='store',
                        type=str, required=False,
                        help="""Output directory to store the results.
                        If not given it assumes automatically the current working
                        directory.""")
    parser.add_argument('-c', '--cutoff', dest='cutoff', action='store',
                        type=float, required=False, default=0.,
                        help="""Cutoff over posteriors. Only values over the cutoff
                        will be added to the sitecount. Default value is 0.""")
    parser.add_argument('-p', '--proxy', dest='proxyBED', action='store',
                        type=str, required=False,
                        help="""Optional BED file that contains regions of interest.
                        Only sites within the regions of interest will be added to the
                        final sitecount. Note that, if not provided, all regions will
                        be considered to build the sitecounts.""")
    args = parser.parse_args()
    return args


def load_sitecounts(infile, cutoff=0., proxy=None):
    """
    :param infile: name of the MotEvo output file for a given motif
    :param cutoff: optional cutoff for the posterior of binding site
    :param proxy: optional BED file containing regions of interest, region ID is at 5th column
    :return: a dictionary for the regions and the position and posterior of sites for each region
    """
    sitecounts = OrderedDict()  # used to remember the insertion order
    if not proxy:    # if the proxy BED file is not provided
        with open(infile, 'r') as inf:
            for rec in csv.reader(inf, delimiter='\t'):
                region_name = rec[3].split(";")[-1]
                post = float(rec[4])
                if post > cutoff:
                    sitecounts.setdefault(region_name, [])
                    chrom = rec[0]
                    start = int(rec[1])
                    end = int(rec[2])
                    strand = rec[10]
                    sitecounts[region_name].append({'post':post,
                                                    'start':start,
                                                    'end':end,
                                                    'chr':chrom,
                                                    'strand':strand})
    else:
        regions_of_interest = dict([(l.split()[4], 0) for l in open(proxy)])   # by default it' assumed that 5th
                                                                               # column contains the region IDs
        with open(infile, 'r') as inf:
            for rec in csv.reader(inf, delimiter='\t'):
                region_name = rec[3].split(";")[-1]
                if not (region_name in regions_of_interest):
                    continue
                post = float(rec[4])
                if post > cutoff:
                    sitecounts.setdefault(region_name, [])
                    chrom = rec[0]
                    start = int(rec[1])
                    end = int(rec[2])
                    strand = rec[10]
                    sitecounts[region_name].append({'post':post,
                                                    'start':start,
                                                    'end':end,
                                                    'chr':chrom,
                                                    'strand':strand})
    return sitecounts


def load_regions_of_interest(fname):
    return dict([(l.split()[4], 0) for l in open(fname)])   # by default it' assumed that 5th column contains the region IDs


def count_nonoverlapping_motifs(sites, filename, cutoff=0.):
    """
    takes the sitecounts for a motif, in a dictionary format, produced by load_sitecount() function,
    and a filename of a MotEvo output file. Then it counts motifs pairs that are overlapping each other.
    :param sites: dictionary type for a given motif
    :param filename: motevo output filename
    ;:param overlap_limit: an integer that indicates how far the two sites need to be from each other
    :param cutoff: minimum cutoff over the posterior
    :param proxy: optional dictionary that contains as key the IDs of regions of interest
    :return: a dictionary that counts double appearance of motifs, where the overlapping sites are filtered
    """
    double_sitecounts = OrderedDict()
    dist = np.zeros(500, dtype=np.float16)
    for region in sites.keys():
        double_sitecounts.setdefault(region, 0)    
    with open(filename, 'r') as inf:
        for rec in csv.reader(inf, delimiter='\t'):
            region_name = rec[3].split(";")[-1]
            double_sitecounts.setdefault(region_name, 0)
            if not (region_name in sites):
                continue
            post = float(rec[4])
            if post > cutoff:
                start = int(rec[1])
                end = int(rec[2])
                for a_site in sites[region_name]:
                    if start >= a_site["start"]:
                        if (a_site["end"] - start) < 0:
                            double_sitecounts[region_name] = 1 
                            distance = int(np.abs((start + (end - start)/2) - \
                                                (a_site["start"] + (a_site["end"] - a_site["start"])/2)))
                            dist[distance] += 1                           
                    else:
                        if (end - a_site["start"]) < 0:
                            double_sitecounts[region_name] = 1
                            distance = int(np.abs((start + (end - start)/2) - \
                                                (a_site["start"] + (a_site["end"] - a_site["start"])/2)))
                            dist[distance] += 1                            
    return double_sitecounts, dist


def main():
    args = arguments()
    motif_sites = load_sitecounts(args.input_file, args.cutoff, args.proxyBED)
    if args.output_dir:
        outdir = args.output_dir
    else:
        outdir = ''
    double_sitecounts = OrderedDict()
    for motif in os.listdir(args.dirname):
        motevo_output_file = os.path.join(args.dirname, motif)
        double_sitecounts.setdefault(motif, 0)
        double_sitecounts[motif], dist = \
                count_nonoverlapping_motifs(motif_sites, motevo_output_file, args.cutoff)
        fname = os.path.join(outdir, '%s.dist' % motif)
        with open(fname, 'w') as outf:
            for i, d in enumerate(dist):
                outf.write('\t'.join([
                    str(i),
                    str(d) + '\n',
                ]))

    with open( os.path.join(outdir, os.path.basename(args.input_file)), 'w') as outf:
        outf.write('\t'.join([motif for motif in double_sitecounts.keys()]) + '\n')
        for region in motif_sites.keys():
            outf.write(region + '\t')
            outf.write('\t'.join([str(double_sitecounts[motif][region]) for motif in double_sitecounts.keys()]))            
            outf.write('\n')
    return 0


if __name__ == '__main__':
    main()

