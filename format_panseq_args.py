import argparse
import os
from subprocess import check_output

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('--fragment', required = True)
    parser.add_argument('--output', required = True)
    parser.add_argument('--querydir', required = True)
    parser.add_argument('--basedir', required = True)
    parser.add_argument('--cores', required = True)
    parser.add_argument('--percent_id', required = True)
    parser.add_argument('--carriage', type = float, required = True)

    return parser.parse_args()

def format_text(querydir, basedir, cores, frag, pid, carriage_thresh):

    def prog_dir(x):
        
        path = check_output(['which', x])
        return os.path.dirname(path) + '/'

    blastdir = prog_dir('blastn') 
    mummerdir = prog_dir('mummer')
    muscle = prog_dir('muscle') + 'muscle'

    args = [
            ('queryDirectory', querydir),
            ('baseDirectory', basedir),
            ('numberOfCores', cores),
            ('mummerDirectory', mummerdir),
            ('blastDirectory', blastdir),
            ('muscleExecutable', muscle),
            ('fragmentationSize', frag),
            ('minimumNovelRegionSize', frag),
            ('percentIdentityCutoff', pid),
            ('coreGenomeThreshold', carriage_thresh),
            ('novelRegionFinderMode', 'no_duplicates'),
            ('maxNumberResultsInMemory', '1000'),
            ('runMode', 'pan')
           ]

    formatted_args = '\n'.join(['\t'.join(x) for x in args])

    return formatted_args


def main():

    args = arguments()

    carriage_thresh = str(int(len(os.listdir(args.querydir)) * (args.carriage / 100)))

    text = format_text(args.querydir,
                       args.basedir,
                       args.cores,
                       args.fragment,
                       args.percent_id,
                       carriage_thresh)

    with open(args.output, 'w') as f:
        f.write(text)

if __name__ == '__main__':
    main()
