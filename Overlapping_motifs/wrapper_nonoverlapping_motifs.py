from os import listdir, path, system
from itertools import combinations

def arguments():
    """adding arguments for the script"""
    import argparse
    parser = argparse.ArgumentParser(description="""
    """)
    parser.add_argument('-d', '--dir', dest='dirname', action='store',
                        type=str, required=True,
                        help='The directory name for sitecounts')
    parser.add_argument('-o', '--dest', dest='destination_dir', action='store',
                        type=str, required=True,
                        help="""The destination directory for output files""")
    args = parser.parse_args()
    return args


def submitJob(input_dir, motif_file, dest):
    """For each file, it submites a job the the queue
    The job type is set to be Normal which means it must take at most
    12 hours to be finished."""
    motif_name = path.basename(motif_file)
    job_name = "%s_nonoverl" % (motif_name)
    dest = path.join(dest, motif_name)
    system('mkdir "%s"' % dest)
    cmd = " ".join([
        "python",
        "~/codes/Motif.Pairs/Overlapping_motifs/count_nonoverlapping_motifs.py \\\n",
        '-d "%s" \\\n' % input_dir,
        '-i "%s" \\\n' % motif_file,
        '-c 0.5 \\\n', 
        '-o "%s" \n' % dest,
    ])    
    with open("job_%s_nonoverlapping.sh" % job_name, "w") as inf:
        inf.write("\n".join([
            '#!/bin/bash',
            '#BSUB -L /bin/bash',
            '#BSUB -o "%s.stdout"' % job_name,
            '#BSUB -e "%s.stderr"' % job_name,            
            '#BSUB -J "%s"' % job_name,
            '#BSUB -M 100000000',
            '#BSUB -R rusage[mem=100000]',            
            '#BSUB normal',
            '',
            cmd,
        ]))
    system('bsub < "%s"' % ("job_%s_nonoverlapping.sh" % job_name))
    return 0


def main():
    args = arguments()
    files = [path.join(args.dirname, f) for f in listdir(args.dirname)]  # temporary for test
    for motif in files:
        dirname = path.join(args.destination_dir, path.basename(motif))
        try:
            resulted_files = listdir(dirname)
        except OSError:
            submitJob(args.dirname, motif, args.destination_dir)
        if not len(resulted_files) > 0:
            submitJob(args.dirname, motif, args.destination_dir)
            print motif
    return 0


if __name__ == '__main__':
    main()
