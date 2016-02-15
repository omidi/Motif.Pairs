from collections import OrderedDict
import csv
from numpy import round

def arguments():
    """adding arguments for the script"""
    import argparse
    parser = argparse.ArgumentParser(description="""
    Takes as input the sitecounts for two motifs,
    as it's produced by "process_motevo_output.py",
    and as output produces a column matrix that each
    element shows the sum of posteriors for motif A
    and motif B in a given region.
    """)
    parser.add_argument('-a', '--motifA', dest='motif_a', action='store',
                        type=str, required=True,
                        help='The first motif file (Note that the order is irrelevant)')
    parser.add_argument('-b', '--motifB', dest='motif_b', action='store',
                        type=str, required=True,
                        help="""The second motif file""")
    parser.add_argument('-o', '--output', dest='output_file', action='store',
                        type=str, required=False,
                        help="""Optional output file name. If not given,
                        the results will be printed in the console.""")
    args = parser.parse_args()
    return args


def count_double_motifs(file_a, file_b):
    """given two motif sitecount files, it extracts
    motif pairs counts"""
    counts1 = dict()
    counts2 = dict()
    regions = OrderedDict()
    with open(file_a) as inf1, open(file_b) as inf2:
        for rec in csv.reader(inf1, delimiter='\t'):
            counts1.setdefault(rec[0], float(rec[1]))
            regions.setdefault(rec[0], 0.)
        for rec in csv.reader(inf2, delimiter='\t'):
            counts2.setdefault(rec[0], float(rec[1]))
            regions.setdefault(rec[0], 0.)

    for region in regions.keys():
        if region in counts1 and region in counts2:
            regions[region] = counts1[region] + counts2[region]

    return regions

def main():
    args = arguments()
    counts = count_double_motifs(args.motif_a, args.motif_b)
    if not args.output_file:
        for region, count in counts.items():
            print '\t'.join([
                region,
                str(round(count, 5)),
            ])
    else:
        with open(args.output_file, 'w') as outf:
            for region, count in counts.items():
                outf.write('\t'.join([
                    region,
                    str(round(count, 5)) + '\n',
                ])
                )

    return 0


if __name__ == '__main__':
    main()
