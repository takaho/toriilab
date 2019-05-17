import os, sys, re
import pandas as pd
import numpy as np

qthr = 0.05
sig = []
with open('THS/quantile_impulse.tsv') as fi:
    fi.readline()
    for line in fi:
        items = line.split('\t')
        # print(items[0])
        if items[2] == 'NA': continue
        pval = float(items[2])
        qval = float(items[3])
        if qval < qthr:
            sig.append(items[0])


#        print(items[0], items[1], pval, qval)

def quantile_normalization(matrix):
    ranks = matrix.stack().groupby(matrix.rank(method='first').stack().astype(int)).mean()
    return matrix.rank(method='min').stack().astype(int).map(ranks).unstack()

def load_normalized(filename, genes=None):
    import hashlib, pickle, gzip
    md5 = hashlib.md5()
    md5.update(os.path.abspath(filename).encode('utf-8'))
    if genes is None:
        md5.update('__all__'.encode('utf-8'))
    else:
        for g in sorted(genes):
            md5.update('_{}'.format(g).encode('utf-8'))
    cache = os.path.join(os.path.dirname(filename), '.' + md5.hexdigest()[0:10] + '.cache')
    if os.path.exists(cache) and os.path.getsize(cache) > 1000:
        with gzip.open(cache, 'rb') as fi:
            return pickle.load(fi)
    # data = {}
    data = quantile_normalization(pd.read_csv(filename, sep='\t', index_col=0))
    if genes is not None:
        data = data.loc[genes]
    with gzip.open(cache, 'wb') as fo:
        pickle.dump(data, fo)
    return data

data = load_normalized('THS/500bp.tsv', sig)

import sklearn.cluster
n_clusters = 6
clu = sklearn.cluster.AgglomerativeClustering(n_clusters=n_clusters)

zscores = data.copy()
nrow, ncol = zscores.shape

zscores = np.divide(np.subtract(zscores, np.mean(data.values, axis=1).reshape(nrow, 1)), np.std(data.values, axis=1).reshape(nrow,1))

print(zscores.shape, np.mean(data, axis=1).shape)
print(zscores)
#zscores /= np.std(data, axis=1)

clu.fit_predict(zscores)

corder = sorted(range(len(sig)), key=lambda i_:clu.labels_[i_])
z = []
for i in corder:
    g = sig[i]
# for g in sig:
    try:
        # values = np.array(data.loc[g])
        # mean = np.mean(values)
        # sd = np.std(values)
        # zscore = (values - mean) / sd
        z.append(zscores.loc[g].values)
    except Exception as e:
        sys.stderr.write(str(e) + '\n')
        continue

import plotly.offline
import plotly.graph_objs as go

hm = go.Heatmap(z=z, y=sig, x=data.columns)
fig = go.Figure([hm, ])
plotly.offline.plot(fig, filename='out.html')
