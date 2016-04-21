import csv
from collections import OrderedDict
import numpy as np
from numpy import log, exp
import os
from scipy.special import gammaln
import re
import random
from string import upper


def arguments():
    """adding arguments for the script"""
    import argparse
    parser = argparse.ArgumentParser(description="""
    Given a set of genes, it generates a FASTA file containing promoters for the given genes
    and then run MotEvo over these DNA sequences using 190 SwissRegulon WMs.
    The output will be located in a directory that user gives at the time of running.
    """)
    parser.add_argument('-i', '--input', dest='genes', action='store',
                        type=str, required=True,
                        help="""A filename that contains at each line name of genes""")
    parser.add_argument('-o', '--output', dest='outputdir', action='store',
                        type=str, required=True,
                        help="""The output directory. If it does not exist before, it will
                        be created.""")
    args = parser.parse_args()
    return args


def main():
    args = arguments()

    os.system('mkdir ' + args.outputdir)
    genes = dict([(x.split()[0].rstrip(), 1) for x in open(args.genes)])

    promoter_annot = '/home/somidi/scratch/Tissue.Specificity/data/promoters/MARA/mara_promoters_gene_name_association.mm10.bed'
    promoter_fasta = '/home/somidi/scratch/Tissue.Specificity/data/promoters/MARA/mm10_mara_promoters.fasta'
    wm_directory = '/home/somidi/scratch/WMs/SwissRegulon'

    promoters = dict()
    with open(promoter_annot) as inf:
        for rec in csv.reader(inf, delimiter='\t'):
            if rec[6] in genes:
                promoters.setdefault(rec[3], rec[6])


    fasta_file = os.path.join(args.outputdir, 'promoters.fa')
    output_promoters = open(fasta_file, 'w')
    with open(promoter_fasta) as inf:
        while True:
            seq_id = inf.readline().rstrip()
            if not seq_id:
                break
            if seq_id.replace('>>', '') in promoters:
                seq = inf.readline().rstrip()
                output_promoters.write('\n'.join([
                    seq_id,
                    upper(seq) + '\n',
                ]))
    output_promoters.close()

    ## now run MotEvo over the newly made promoters set

    cmd = ' '.join([
        'python',
        'calculate_sitecount.py',
        '-w %s' % wm_directory,
        '-f %s' % fasta_file,
        '-o %s' % args.outputdir,
    ])
    os.system(cmd)





if __name__ == '__main__':
    main()