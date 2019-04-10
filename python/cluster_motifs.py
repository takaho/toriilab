"""
Clustering motifs
"""
import numpy as np
import pandas as pd

import sklearn.cluster
import plotly.offline
import plotly.graph_objs as go
import sys, argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i' )
parser.add_argument('-n', default=7, type=int )
parser.add_argument('-o', default=None)
args = parser.parse_args()

motifs = pd.read_csv(args.i, sep='\t', index_col=0)
tfs = motifs.columns[1:]
matrix = motifs[tfs].copy()
matrix[matrix < 1] = 1
n_clusters = args.n
clu = sklearn.cluster.AgglomerativeClustering(n_clusters=n_clusters)
clu.fit_predict(np.log(matrix))
# print(len(clu.labels_))

clorder = [[] for i in range(n_clusters)]
values = matrix.values
ethr = 2
cluster_score = [[] for i in range(n_clusters)]
n_cols = matrix.shape[1]
for cl in range(n_clusters):
    scores = [0] * n_cols
    n = 0
    for i, label in enumerate(clu.labels_):
        if label != cl: continue
        for j, v in enumerate(matrix.values[i]):
            if v > ethr:
                scores[j] += 1
        n += 1
    amx = np.argmax(scores)
    if n == 0: n += 1
    if scores[amx] < n / 4:
        score = 10
    else:
        score = amx
    clorder[cl] = score

motif_names= motifs.index
tf_names = matrix.columns
ost = sys.stdout
if args.o is not None:
    ost = open(args.o, 'w')
ost.write('Motif\tCluster\tSequence\t' + '\t'.join(tfs) + '\n')
for cl in sorted(range(n_clusters), key=lambda v:clorder[v]):
    for i, label in enumerate(clu.labels_):
        if label != cl:
            continue
        name = motif_names[i]
        ostr = '{}\t{}\t{}'.format(name, label, motifs['Motif'][name])
        vmax = 0
        for v in matrix.loc[name].values:
            if v <= 1.0:
                ostr += '\t.'
            else:
                vmax = max(vmax, v)
                ostr += '\t{:.2f}'.format(np.log2(v))
        if 1 or vmax > ethr: ost.write(ostr + '\n')
if args.o: ost.close()