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
parser.add_argument('-g', default='out.pdf')
args = parser.parse_args()

motifs = pd.read_csv(args.i, sep='\t', index_col=0)
tfs = motifs.columns[1:]
matrix = motifs[tfs].copy()
matrix[matrix < 1] = 1
n_clusters = args.n
clu = sklearn.cluster.AgglomerativeClustering(n_clusters=n_clusters)
clu.fit_predict(np.log(matrix))

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

if args.g:
    import reportlab.pdfgen.canvas
    import reportlab.lib.colors
    cnv = reportlab.pdfgen.canvas.Canvas(args.g)
    cnv.setFont('Helvetica', 8)
    left, top = 200, 700
    xstep, ystep = 25, 8
    # SPCH 4,205,4
    # MUTE 3,205,255
    # FAMA 154,2,255

    colors = [
        reportlab.lib.colors.Color(.8, .1, .3, .5), # AtML1
        reportlab.lib.colors.Color(.02, 1., .02, .5), # SPCH
        reportlab.lib.colors.Color(.02, .8, 1., .5), # MUTE
        reportlab.lib.colors.Color(.6, .01, 1., .5), # FAMA
        reportlab.lib.colors.Color(.3, .6, .1, .5), # GC1
    ]
    gray = reportlab.lib.colors.Color(0,0,0,.4)
    for i, tf in enumerate(tfs):
        cnv.setFillColor('black')
        cnv.drawString(left + xstep * (i + .5) - cnv.stringWidth(tf) * .5, top + 15, tf)
else:
    cnv = None

y = top
for cl in sorted(range(n_clusters), key=lambda v:clorder[v]):
    for i, label in enumerate(clu.labels_):
        if label != cl:
            continue
        name = motif_names[i]
        ostr = '{}\t{}\t{}'.format(name, label, motifs['Motif'][name])
        if cnv:
            cnv.setFillColor('black')
            cnv.drawString(left - 70, y, name)
            cnv.drawString(left + xstep * (len(tfs) + .5), y, motifs['Motif'][name])
        vmax = 0
        for j, enrichment in enumerate(matrix.loc[name].values):
            if enrichment <= 1.0:
                ostr += '\t.'
            else:
                vmax = max(vmax, enrichment)
                ostr += '\t{:.2f}'.format(np.log2(enrichment))
            if cnv:
                cx = left + xstep * (j + .5)
                cy = y + ystep * .5
                xsize = (.2 + min(5, np.log2(enrichment))) * xstep * .15
                ysize = xsize / xstep * ystep * 2
                if enrichment <= 1.0:
                    cnv.setFillColor(gray)
                    roundlevel = 0.5
                else:
                    cnv.setFillColor(colors[j])
                    roundlevel = xsize / 2
                cnv.roundRect(cx - xsize, cy - ysize, xsize * 2, ysize * 2, roundlevel, 0, fill=1)
                    #cnv.rect(cx + xstep * j - xsize, cy - ysize, xsize * 2, ysize * 2, 0, fill=1)
        y -= ystep
        if 1 or vmax > ethr: ost.write(ostr + '\n')
if args.o: ost.close()
if cnv:
    for i in range(0, 5):
        logenr = i * .8 + 1
        xsize = (.2 + min(5, logenr)) * xstep * .15
        ysize = xsize / xstep * ystep * 2
        roundlevel = xsize / 2
        cnv.setFillColor('navy')
        cnv.roundRect((i  * 2)  * xstep * .5 + left - xsize, y - ysize - 100, xsize * 2, ysize * 2, roundlevel, 0, 1)
    cnv.save()