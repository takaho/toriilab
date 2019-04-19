import hashlib, gzip, pickle
import os, sys, re

def load_peaks(filename, upstream=500, downstream=0, nuconly=False, cache=None):
    if cache is None:
        md5 = hashlib.md5()
        md5.update(os.path.abspath(filename).encode('utf-8'))
        md5.update('::{}::{}::{}'.format(upstream, downstream, int(nuconly)).encode('utf-8'))
        cache = os.path.join(os.path.dirname(filename),  '.' + md5.hexdigest() + '.cache')
    if os.path.exists(cache) and os.path.getsize(cache) > 10000:
        with gzip.open(cache) as fz:
            return pickle.load(fz)
    data = set()
    with open(filename) as fi:
        for line in fi:
            items = line.strip().split('\t')
            for elem in items[2].split('///'):
                elems = re.split('\\s*//\\s*', elem.strip())
                if len(elems) >= 4:
                    loc, dist = elems[3].split(':')
                    dist = int(dist)
                    if (loc == 'UP' and dist <= upstream) or (loc == 'IN' and dist <= downstream):
                        # print(loc, dist)
                        data.add(elems[0])
    with gzip.open(cache, 'wb') as fz:
        pickle.dump(data, fz)
    return data

def load_gene_location(gtf, cache=None):
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

import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-i', nargs='+' )
parser.add_argument('--upstream', type=int, default=500)
parser.add_argument('--downstream', type=int, default=0)
parser.add_argument('--gtf', default='/mnt/smb/tae/tair10/Arabidopsis_thaliana.TAIR10.40.chrm.gtf')
parser.add_argument('--verbose', action='store_true')
parser.add_argument('-n', type=int, default=10)
args = parser.parse_args()
genes = load_gene_location(args.gtf)
# print(list(genes.values())[0:10])
num_total = len([g for g in genes.values() if re.match('chr\\d', g[1])])
# print(num_total)
samples = []
data = {}
used_genes = set()
for fn in args.i:
    basename = os.path.basename(fn)
    if basename.find('.') >= 0:
        name = basename[:basename.rfind('.')]
    else:
        name = basename
    genes = load_peaks(fn, args.upstream, args.downstream, True)
    print(name, len(genes))
    if len(genes) > 0:
        samples.append(name)
        data[name] = genes
        for g in genes: used_genes.add(g)
N = len(used_genes)
M = len(samples)
matrix = np.zeros((N, M), dtype=np.float)
genes = list(sorted(used_genes))
vals = []
for i, s in enumerate(samples):
    datum = data[s]
    for j, g in enumerate(genes):
        vals.append(int(g in datum))
df = pd.DataFrame(np.array(vals).reshape(N, M), index=genes, columns=samples)

for i, si in enumerate(samples):
    for j in range(i + 1, len(samples)):
        sj = samples[j]
        n_0 = len(data[si])
        n0_ = len(data[sj])
        n00 = 0
        for g in data[si]:
            if g in data[sj]:
                n00 += 1
        n01 = n_0 - n00
        n10 = n0_ - n00
        n11 = num_total - n00 - n01 - n10
        n_1 = n01 + n11
        n1_ = n10 + n11
        jaccard = n00 / (n_0 + n0_ - n00)
        matrix[i][j] = matrix[j][i] = jaccard
        # print(i, j, n00, n01, n10, n11)
        odds = float(n00 * n11) / (n01 * n10)
        phi = float(n11 * n00 - n01 * n10) / np.sqrt(n_0 * n0_ * n_1 * n1_)
        print(i, j, jaccard, odds, phi)

print(matrix)

import scipy.cluster.hierarchy
Z = scipy.cluster.hierarchy.linkage(df.T, method='ward', metric='euclidean')
print(Z)
exit()    

print(df.shape)
import sklearn.cluster
n_clusters = args.n
ag = sklearn.cluster.AgglomerativeClustering(n_clusters=n_clusters)
ag.fit_predict(df)
results = np.zeros((n_clusters, M))
nums = [0] * n_clusters
matrix = df.values

def load_tpm(cache='.tpm.cache'):
    if os.path.exists(cache) and os.path.getsize(cache) > 1000:
        with gzip.open(cache) as fi:
            return pickle.load(fi)
    data = {}
    with open('SRP043607.count') as fi:
        samples = fi.readline().strip().split('\t')
        for line in fi:
            items = line.strip().split('\t')
            gene = items[0].split('.')[0]
            vals = [float(x) for x in items[1:]]
            if gene in data:
                row = data[gene]
                for i,v  in enumerate(vals):
                    row[i] += v
            else:
                data[gene] = vals
    matrix = []
    indexes = []
    for gene, vals in data.items():
        indexes.append(gene)
        matrix.append([(vals[0] + vals[1]) * .5, (vals[2] + vals[3] + vals[4]) / 3, (vals[5] + vals[6]) * .5, (vals[7] + vals[8]) * .5, (vals[9] + vals[10]) * .5])
    df = pd.DataFrame(matrix, index=indexes, columns=['ML1', 'SPCH', 'MUTE', 'FAMA', 'GC1'])
    with gzip.open(cache, 'wb') as fi:
        pickle.dump(df, fi)
    return df
 
import pandas as pd
tpm = load_tpm()
mintpm = 0.1
table = tpm.copy()
table[table < mintpm] = mintpm
logtpm = np.log2(table)
available_genes = set(table.index)
cl2genes = [[] for i in range(n_clusters)]
expr_columns = table.columns
num_expr = len(expr_columns)
for i in range(N):
    label = ag.labels_[i]
    results[label] += matrix[i]
    nums[label] += 1    
    cl2genes[label].append(genes[i])

for label in range(n_clusters):
    expr = np.zeros(num_expr)
    n = 0
    for g in cl2genes[label]:
        if g in available_genes:
            expr += logtpm.loc[g].values
            n += 1
    estr = '['
    for v in expr:
        estr += ' {:.2f}'.format(v / n)
    estr += ']'
    print('{}\t{}\t{}\t{}'.format(label, nums[label], results[label] / nums[label], estr))


exit()

for i, si in enumerate(samples):
    for j, sj in enumerate(samples):
        if i < j:
            matrix[j][i] = matrix[i][j]
        else:
            n_0 = len(si)
            n0_ = len(sj)
            n00 = 0
            for g in si:
                if g in sj:
                    n00 += 1
            n01 = n_0 - n00
            n10 = n0_ - n00
            n11 = num_total - n00 - n01 - n10
            jaccard = n00 / (n_0 + n0_ - n00)
            matrix[i][j] = jaccard

print(scipy.cluster.hierarchy.dendrogram(Z))



