import os, sys, re
import numpy as np
import argparse

atac_columns = [[3,4,5], [6,7,8], [9,10,11,12,13], [14,15,16], [17,18,19]]
data = []
cluster = []
markers = []
marker_genes = 'ATML1', 'SPCH', 'MUTE', 'FAMA', 'GASA9'
with open(sys.argv[1]) as fi:
    items = fi.readline().strip().split('\t')
    columns = items[3:]
    if len(columns) > 10:
        columns = list(marker_genes) + columns[-5:]
    locations = []
    for line in fi:
        avals = []
        evals = []
        items = line.strip().split('\t')
        for cols in atac_columns:
            avals.append(np.mean([float(items[col]) for col in cols]))
        for item in items[-5:]:
            if item == '.':
                evals = None
                break
            evals.append(np.log(float(item) + .1))
        if avals is not None and evals is not None:
            for m in marker_genes:
                if items[1].find('({})'.format(m)) >= 0:
                    markers.append((len(locations), m))
                    # print(m)
                    break
            locations.append(items[1])
            data.append(avals + evals)
            # print(avals, evals)
            cluster.append(4 - int(items[0]))

import umap
import pandas as pd
import sklearn.decomposition


matrix = pd.DataFrame(data, columns=columns, index=locations)
u = umap.UMAP(n_components=3)
res = u.fit_transform(matrix)

# ipca = sklearn.decomposition.IncrementalPCA(n_components=3)
# res = ipca.fit_transform(matrix)

#print(res.shape)
import plotly
import plotly.graph_objs as go
s = 4
traces = [go.Scatter3d(x=res[:,0], y=res[:,1], z=res[:,2], text=matrix.index, hoverinfo='text', mode='markers', marker=dict(size=s, colorscale='Viridis', cmin=0, cmax=max(cluster), opacity=0.8, color=cluster))]

layout = go.Layout(margin=dict(l=0, r=0, b=0, t=0), scene=dict(xaxis=dict(title='UMAP1'), yaxis=dict(title='UMAP2'), zaxis=dict(title='UMAP3')))

for i, m in markers:
    pos = res[i:i+1].T
    # print(cluster[i:i+1])
    print(i, cluster[i], m)
    trace = go.Scatter3d(x=pos[0], y=pos[1], z=pos[2], mode='markers', marker=dict(size=s * 2, colorscale='Viridis', color=[cluster[i],], cmin=0, cmax=max(cluster), line=dict(width=2, color='rgb(192,64,128)'), symbol='diamond'), name=m)
    traces.append(trace)

fig = go.Figure(data=traces, layout=layout)
plotly.offline.plot(fig, filename='scatter.html')
#for i in range(matrix.shape[0]):

