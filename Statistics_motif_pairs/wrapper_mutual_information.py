from os import listdir, path, system
from itertools import combinations

def arguments():
    """adding arguments for the script"""
    import argparse
    parser = argparse.ArgumentParser(description="""
    """)
    parser.add_argument('-r', '--regions', dest='regions', action='store',
                        type=str, required=True,
                        help='The regions of interest in BED format')
    parser.add_argument('-o', '--dest', dest='destination_dir', action='store',
                        type=str, required=True,
                        help="""The destination directory for output files""")
    parser.add_argument('-c', '--cutoff', dest='cutoff', action='store',
                        type=float, required=False, default=0.5,
                        help="""Minimum posterior for selecting TFBS""")
    args = parser.parse_args()
    if args.cutoff > 1 or args.cutoff<0:
        print 'Provided cutoff is invalid.'
        print 'Program halts!'
        exit()
    return args


def submitJob(motif_file, dest, regions, cutoff=0.5):
    """For each file, it submites a job the the queue
    The job type is set to be Normal which means it must take at most
    12 hours to be finished."""
    motif_name = path.basename(motif_file)
    job_name = "%s_mutInfo" % (motif_name)
    dest = path.join(dest, motif_name)
    cmd = " ".join([
        "python",
        "/home/somidi/scratch/Tissue.Specificity/codes/Motif.Pairs/Statistics_motif_pairs/mutual_information_motif_pairs.py \\\n",
        '-m "%s" \\\n' % motif_name,
        '-r "%s" \\\n' % regions,
        '-c %f > %s' % (cutoff, path.join(dest, motif_name)),
    ])    
    with open("job_%s_mutual_info.sh" % job_name, "w") as inf:
        inf.write("\n".join([
            '#!/bin/bash',
            '#BSUB -L /bin/bash',
            '#BSUB -o "%s.stdout"' % job_name,
            '#BSUB -e "%s.stderr"' % job_name,
            '#BSUB -J "%s"' % job_name,
            # '#BSUB -M 100000000',
            '#BSUB -R rusage[mem=250000]',
            '#BSUB normal',
            '',
            cmd,
        ]))
    system('bsub < "%s"' % ("job_%s_mutual_info.sh" % job_name))
    return 0


def main():
    args = arguments()
    motevo_directory = "/home/somidi/scratch/Tissue.Specificity/Motif.Pairs/MotEvo.Outputs/min.post.50"
    files = [path.join(motevo_directory, f) for f in listdir(motevo_directory)]
    for motif in files:
        submitJob(motif, args.destination_dir, args.regions, args.cutoff)
    return 0


if __name__ == '__main__':
    main()
