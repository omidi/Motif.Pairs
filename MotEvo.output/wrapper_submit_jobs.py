from os import listdir, path, system

def arguments():
    """adding arguments for the script"""
    import argparse
    parser = argparse.ArgumentParser(description="""
    Given a directory that contains all the MotEvo output files.
    NOTE: make sure that the input directory only contains the MotEvo output files.
    """)
    parser.add_argument('-d', '--dir', dest='dirname', action='store',
                        type=str, required=True,
                        help='The directory name for MotEvo output files.')
    parser.add_argument('-o', '--dest', dest='destination_dir', action='store',
                        type=str, required=True,
                        help="""The destination directory to contain the output files""")
    parser.add_argument('-c', '--cutoff', dest='cutoff', action='store',
                        type=float, required=False, default=0.,
                        help="""Optional cutoff over the posteriors for each motif. It
                        will be the same for all motifs in the directory.
                        The default value is set to zero.""")
    args = parser.parse_args()
    return args


def submitJob(filename, dest, cutoff=0.):
    """For each file, it submites a job the the queue
    The job type is set to be Normal which means it must take at most
    12 hours to be finished."""
    name = path.basename(filename)
    cmd = " ".join([
        "process_motevo_output.py",
        "-i %s" % filename,
        "-o %s" % path.join(dest, name),
        "-c %f" % cutoff,
    ])
    with open("job_%s.sh" % name) as inf:
        inf.write("\n".join([
            '#!/bin/bash',
            '#BSUB -L /bin/bash',
            '#BSUB -o "%s.stdout"' % name,
            '#BSUB -e "%s.stderr"' % name,
            '#BSUB -J "%s"' % name,
            '#BSUB normal',
            '',
            cmd,
        ]))
    # system('bsub < %s' % ("job_%s.sh" % name))
    return 0


def main():
    args = arguments()
    files = [path.join(args.dirname, f) for f in listdir(args.dirname)]
    for file in files:
        submitJob(file, args.destination_dir, args.cutoff)
    return 0


if __name__ == '__main__':
    main()