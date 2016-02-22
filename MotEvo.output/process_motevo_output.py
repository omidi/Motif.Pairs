import csv
from collections import OrderedDict
from numpy import round

def arguments():
    """adding arguments for the script"""
    import argparse
    parser = argparse.ArgumentParser(description="""
    takes an output file from MotEvo for a given motif,
    and produces a column matrix where each element is the
    sum of posteriors for a given region.
    """)
    parser.add_argument('-i', '--input', dest='input_file', action='store',
                        type=str, required=True,
                        help='The MotEvo output file')
    parser.add_argument('-output', '--output', dest='output_file', action='store',
                        type=str, required=False,
                        help="""Optional output file name, if not given the output
                        will be printed out in the console""")
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


def makeSitecount(infile, outfile=None, cutoff=0., proxy=None):
    """it takes the name of a MotEvo file and returns a dictionary with
    that contains the sum of posterior for each region
    """
    sitecounts = OrderedDict()  # used to remember the insertion order
    if not proxy:    # if the proxy BED file is not provided
        with open(infile, 'r') as inf:
            for rec in csv.reader(inf, delimiter='\t'):
                region_name = rec[3].split(";")[-1]
                sitecounts.setdefault(region_name, 0.)
                post = float(rec[4])
                if post > cutoff:
                    sitecounts[region_name] += post
    else:
        regions_of_interest = dict([(l.split()[4], 0) for l in open(proxy)])   # by default it' assumed that 5th
                                                                               # column contains the region IDs
        with open(infile, 'r') as inf:
            for rec in csv.reader(inf, delimiter='\t'):
                region_name = rec[3].split(";")[-1]
                if not (region_name in regions_of_interest):
                    continue
                sitecounts.setdefault(region_name, 0.)
                post = float(rec[4])
                if post > cutoff:
                    sitecounts[region_name] += post
    return sitecounts


def main():
    args = arguments()
    sitecounts = makeSitecount(args.input_file, args.output_file, args.cutoff, args.proxyBED)
    if not args.output_file:
        for region, sitecount in sitecounts.items():
            print '\t'.join([
                region,
                str(round(sitecount, 5)),
            ])
    else:
        with open(args.output_file, 'w') as outf:
            for region, sitecount in sitecounts.items():
                outf.write("\t".join([
                region,
                str(round(sitecount, 5)) + '\n',
                ]))
    return 0


if __name__ == '__main__':
    main()


