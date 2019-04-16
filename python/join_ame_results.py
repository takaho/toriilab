import os, sys, re

class AMEResult(object):
    def __init__(self, motif_id, consensus, pvalue, qvalue, evalue, hits, positive_ratio, negative_ratio):
        self.motif_id = motif_id
        self.consensus = consensus
        self.pvalue = pvalue
        self.qvalue = qvalue
        self.evalue = evalue
        self.hits = hits
        self.positive = positive_ratio
        self.negative = negative_ratio

    def __repr__(self):
        return('{}\t{}\t{:.0f}\t{:.3f}'.format(self.motif_id, self.consensus, self.hits, self.enrichment))
    enrichment = property(lambda s:s.positive / s.negative)
    @classmethod
    def load_ame(cls, filename, qvalue_threshold=1.001):
        results = []
        with open(filename) as fi:
            fi.readline()
            for line in fi:
                items = line.strip().split('\t')
                if len(items) >= 12:
                    qvalue = float(items[7])
                    if qvalue_threshold > qvalue:
                        res = AMEResult(items[2], items[4], float(items[5]), float(items[7]), float(items[6]), int(items[13]), float(items[14]), float(items[16]))
                        # print(repr(res))
                        results.append(res)
        return results
import collections
import openpyxl
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', nargs='+', help='input AME files or their directories')
parser.add_argument('-n', type=int, default=8, help='clusters')
parser.add_argument('--log', action='store_true')
parser.add_argument('--verbose', action='store_true')
parser.add_argument('--evalue', type=float, default=0, help='e-value threshold')
parser.add_argument('--qvalue', type=float, default=1.0, help='q-value threshold')
args = parser.parse_args()
verbose = args.verbose
results = collections.OrderedDict()
motifs = {}
for fn in args.i:
    basename = os.path.basename(fn)
    if basename == 'ame.tsv':
        filename = fn
        name = os.path.basename(os.path.dirname(fn)).split('.')[0]
    else:
        name = os.path.basename(fn).split('.')[0]
        filename = os.path.join(fn, 'ame.tsv')
    res = AMEResult.load_ame(filename, qvalue_threshold=2)
    accepted = {}
    for r in res:
        if (args.qvalue >= 1.0 or r.qvalue < args.qvalue) and (args.evalue <= 0 or args.evalue < r.evalue):
            motifs[r.motif_id] = r.consensus
            if r.motif_id not in accepted:
                accepted[r.motif_id] = r
    results[name] = accepted
    if verbose:
        sys.stderr.write('Loaded {}\t{}\t{}\n'.format(name, len(res), len(accepted)))

import numpy as np
import pandas as pd
# motifs = list(sorted(motifs))
motif_ids = list(sorted(motifs.keys()))
sys.stderr.write('motifs :{}\n'.format(len(motif_ids)))
N = len(motifs)
M = len(results)

# matrix = np.full((N, M), 0.0)
matrix = []
import sklearn.cluster
for i, m in enumerate(motif_ids):
    row = [0] * M
    for j, n in enumerate(results.keys()):
        if m in results[n]:
            r = results[n][m]
        # for r in results[n]:
        #     if r.motif_id == m:
                # print(r.enrichment)
            # print(r.motif_id, m)
            if args.log:
                row[j]= np.log2(r.enrichment)
            else:
                row[j]= r.enrichment - 1#np.log2(r.enrichment)
                # break
    matrix.append(row)
    # print(matrix[i])
df = pd.DataFrame(matrix, index=motif_ids, columns=results.keys())
# print(df)
# print(matrix)
# exit()
n_clusters = args.n
clu = sklearn.cluster.AgglomerativeClustering(n_clusters=n_clusters)#, linkage='complete'
# clu = sklearn.cluster.KMeans(n_clusters=n_clusters, algorithm='elkan')
# .AgglomerativeClustering(n_clusters=n_clusters
# clu = sklearn.cluster.KMeans(n_clusters=n_clusters)
clu.fit_predict(matrix)
# print(len(clu.labels_))

ordered_labels = list(range(n_clusters))
header = 'Cluster\tMotif\tConsensus\t' + '\t'.join(results.keys())
print(header)
stats = []
nums = []
for cluster in ordered_labels:
    aggr = [[] for i in range(M)]
    # averages = [0] * M
    n = 0
    for i, m in enumerate(motif_ids):
        if clu.labels_[i] == cluster:
            ostr = '{}\t{}\t{}'.format(cluster, m, motifs[m])
            for j, enr in enumerate(matrix[i]):
                aggr[j].append(enr)
                # averages[j] += enr
                ostr += '\t{:.2f}'.format(enr)
                n += 1
            print(ostr)
    stats.append([np.median(x_) for x_ in aggr])
    # stats.append([x_ / max(n, 1) for x_ in averages])
    nums.append(n)
for i, avg in enumerate(stats):
    ost = '{}\t{}'.format(i, nums[i])
    for v in avg:
        ost += '\t{:.3f}'.format(v)
    print(ost)
    

