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
                        The default value is set to 0.5""")
    parser.add_argument('-m', '--overlap', dest='overlap', action='store',
                        type=int, required=False, default=0,
                        help="""Determines the limit within which two sites considered
                        overlapping. By default is set to zero, implying two sites are side
                        by side by zero base""")        
    parser.add_argument('-p', '--proxy', dest='proxyBED', action='store',
                        type=str, required=False,
                        help="""A BED file that its 5th column contains the ID of
                        regions of interest. If not given, sitecounts for all the
                        regions within the MotEvo output file will be given""")    
    args = parser.parse_args()
    return args


def submitJob(filename, dirname, dest, cutoff=0.5, proxy=None, overlap=0):
    """For each file, it submites a job the the queue
    The job type is set to be Normal which means it must take at most
    12 hours to be finished."""
    name = path.basename(filename)
    dest_dir = path.join(dest, name)
    cmd = 'mkdir "%s"\n' % dest_dir 
    cmd += " ".join([
        "python",
        "~/codes/Motif.Pairs/MotEvo.output/one_motif_vs_all.py",
        '-i "%s" \\\n' % filename,
        '-d "%s" \\\n' % dirname,
        '-o "%s" \\\n' % dest_dir,
        '-p "%s" \\\n' % proxy if proxy else '',
        '-m %d \\\n' % overlap,         
        '-c %f \n' % cutoff,
    ])
    # cmd += "cd %s\n" % dest 
    # cmd += 'gzip -f "%s"\n' % name
    
    with open("job_%s.sh" % name, "w") as inf:
        inf.write("\n".join([
            '#!/bin/bash',
            '#BSUB -L /bin/bash',
            '#BSUB -o "%s.stdout"' % name,
            '#BSUB -e "%s.stderr"' % name,
            '#BSUB -J "%s"' % name,
            '#BSUB -M 1000000',
            '#BSUB -R rusage[mem=1000]',
            '#BSUB normal',
            '',
            cmd,
        ]))
    system('bsub < "%s"' % ("job_%s.sh" % name))
    return 0


def main():
    args = arguments()
    files = [path.join(args.dirname, f) for f in listdir(args.dirname)]
    for file in files:
        submitJob(file, args.dirname, args.destination_dir, args.cutoff, args.proxyBED, args.overlap)
    return 0


if __name__ == '__main__':
    main()
