import os, sys, re
import argparse

def load_genome(fn):
    genome = {}
    with open(fn) as fi:
        chrom = seq = None
        for line in fi:
            if line.startswith('>'):
                if chrom is not None:
                    genome[chrom] = seq
                chrom = line[1:].strip()
                seq = ''
            elif chrom is not None:
                seq += line.strip()
        if chrom is not None:
            genome[chrom] = seq
    return genome

parser = argparse.ArgumentParser()
parser.add_argument('-i', nargs='+')
parser.add_argument('-g', default='/mnt/smb/tae/tair10/tair10.fa')
parser.add_argument('-s', type=int, default=500)

args = parser.parse_args()
genome = None

if args.i.endswith('.gtf'):

elif args.i.endswith('.bed') or args.i.endswith('.xls'):
    size = args.s
    for fn in args.i:
        name = os.path.basename(fn).split('.')[0]
        sys.stderr.write(name + '\n')
        fn_out = name + '.fa'
        if genome is None:
            genome = load_genome(args.g)
        with open(fn) as fi, open(fn_out, 'w') as fo:
            for line in fi:
                items = line.strip().split('\t')
                # items = line.decode('ascii').strip().split('\t')
                seq = genome.get(items[0], '')
                center = min(len(seq) - size // 2, max(size // 2, (int(items[1]) + int(items[2])) // 2))
                start = center - size // 2
                stop = start + size
                frag = seq[start:stop]
                if len(frag) > size // 2:
                    fo.write('>{}_{}:{}-{}\n{}\n'.format(name, items[0], start, stop, frag))

