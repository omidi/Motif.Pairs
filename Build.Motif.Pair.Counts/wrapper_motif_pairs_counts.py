from os import listdir, path, system
from itertools import combinations

def arguments():
    """adding arguments for the script"""
    import argparse
    parser = argparse.ArgumentParser(description="""
    Given a directory that contains all the motif sitecounts,
    as produced by "process_motevo_output.py" script, this script
    submits for each pairs of motifs a job to the queue to calculate
    the pair-counts (sum of two posteriors).
    """)
    parser.add_argument('-d', '--dir', dest='dirname', action='store',
                        type=str, required=True,
                        help='The directory name for sitecounts')
    parser.add_argument('-o', '--dest', dest='destination_dir', action='store',
                        type=str, required=True,
                        help="""The destination directory for output files""")
    args = parser.parse_args()
    return args


def submitJob(file1, file2, dest):
    """For each file, it submites a job the the queue
    The job type is set to be Normal which means it must take at most
    12 hours to be finished."""
    name1 = path.basename(file1)
    name2 = path.basename(file2)
    job_name = "%s_with_%s" % (name1, name2)
    cmd = " ".join([
        "python",
        "~/codes/Motif.Pairs/Build.Motif.Pair.Counts/motif_pairs_counts.py",
        '-a "%s" \\\n' % file1,
        '-b "%s" \\\n' % file2,
        '-o "%s" \\\n' % path.join(dest, job_name),
    ])
    cmd += "cd %s\n" % dest
    cmd += 'gzip -f "%s"\n' % job_name
    with open("job_%s.sh" % job_name, "w") as inf:
        inf.write("\n".join([
            '#!/bin/bash',
            '#BSUB -L /bin/bash',
            '#BSUB -o "%s.stdout"' % job_name,
            '#BSUB -e "%s.stderr"' % job_name,
            '#BSUB -J "%s"' % job_name,
            '#BSUB -M 500000',
            '#BSUB -R rusage[mem=500]',
            '#BSUB normal',
            '',
            cmd,
        ]))
    system('bsub < "%s"' % ("job_%s.sh" % job_name))
    return 0


def main():
    args = arguments()
    files = [path.join(args.dirname, f) for f in listdir(args.dirname)]
    for pairs in combinations(files, 2):  # for every pair but the order is irrelevant
        submitJob(pairs[0], pairs[2], args.destination_dir)
    return 0


if __name__ == '__main__':
    main()
