import argparse
import os, sys, re, hashlib, gzip, pickle
# import numpy as np
# import pandas as pd

def load_gene_location(gtf, cache=None):
    """ return dict {key:gene_id, value: (name, chrom, start, stop, orientation) }
    """
    if cache is None:
        md5 = hashlib.md5()
        md5.update(os.path.abspath(gtf).encode('utf-8'))
        # md5.update('::{}'.format(int(nuconly).encode('utf-8')))
        cache = os.path.join(os.path.dirname(gtf), '.' + md5.hexdigest() + '.cache')
    if os.path.exists(cache) and os.path.getsize(cache) > 10000:
        with gzip.open(cache) as fz:
            return pickle.load(fz)
    genes = {}
    with open(gtf) as fi:
        for line in fi:
            items = line.strip().split('\t')
            if items[2] != 'gene': continue
            m = re.search('gene_id\\s+"([^"]+)', items[8])
            if m:
                gene_id = m.group(1)
                if gene_id in genes: continue
            else:
                continue
            m = re.search('gene_name\\s+"([^"]+)', items[8])
            if m:
                gene_name = m.group(1)
            else:
                gene_name = gene_id
            genes[gene_id] = gene_name, items[0], int(items[3]), int(items[4]), items[6]
    with gzip.open(cache, 'wb') as fz:
        pickle.dump(genes, fz)
    return genes

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # parser.add_argument('-i', nargs='+' )
    parser.add_argument('--upstream', type=int, default=500)
    parser.add_argument('--downstream', type=int, default=0)
    parser.add_argument('--gtf', default='/mnt/smb/tae/tair10/Arabidopsis_thaliana.TAIR10.40.chrm.gtf')
    parser.add_argument('--verbose', action='store_true')
    # parser.add_argument('-n', type=int, default=10)
    args = parser.parse_args()
    genes = load_gene_location(args.gtf)
    upstream = args.upstream    
    downstream = args.downstream

    for gene, loc in genes.items():
        gene_name, chrom, start, stop, strand = loc
        if gene != gene_name:
            gene_name = '{}({})'.format(gene, gene_name)
        else:
            gene_name = gene
        if strand == '+':
            print('{}\t{}\t{}\t{}'.format(chrom, start - upstream, start + downstream, gene_name))
        else:
            print('{}\t{}\t{}\t{}'.format(chrom, stop - downstream, stop + upstream, gene_name))
