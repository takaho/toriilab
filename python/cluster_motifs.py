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
celltypes = motifs.columns[1:]
matrix = motifs[celltypes].copy()
matrix[matrix < 1] = 1
cmat = matrix.copy()
cmat[cmat <= 1] = 0
cmat[cmat > 2] = 1
n_clusters = args.n
clu = sklearn.cluster.AgglomerativeClustering(n_clusters=n_clusters)
# clu.fit_predict(np.log2(matrix))
clu.fit_predict(cmat)#np.log2(matrix))

# clorder = [[] for i in range(n_clusters)]
values = matrix.values
ethr = 1
cluster_score = [[] for i in range(n_clusters)]
n_cols = matrix.shape[1]
smat = []
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
    # print(cl, scores, amx, n)
    if n == 0: n = 1
    smat.append([cl, ] + [x_ / n for x_ in scores])
index = 0
def __gen_key(a):
    n = len(a) - 1
    value = 0
    for i in range(n):
        value *= 100
        value += min(10, a[i + 1])
        # if a[i] != b[i]:
        #     return -(a[i] < b[i])
    return value
# clorder = [0] * len(smat)
# print(smat)
# for i, l in enumerate(sorted(smat, key=__gen_key)):
#     clorder[l[0]] = i
#     print(__gen_key(l))

# print(clorder)
clorder = [0, 1, 2, 3, 4]
# clorder = 
motif_names= motifs.index
tf_names = matrix.columns
ost = sys.stdout
if args.o is not None:
    ost = open(args.o, 'w')
ost.write('Motif\tCluster\tSequence\t' + '\t'.join(celltypes) + '\n')
print(ost)
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
    for i, celltype in enumerate(celltypes):
        cnv.setFillColor('black')
        cnv.drawString(left + xstep * (i + .5) - cnv.stringWidth(celltype) * .5, top + 15, celltype)
else:
    cnv = None

y = top

def enrichment2size(enr, maxval=8):
    if enr <= 1.0:
        return 0.1
    else:
        return (1.0 + np.log2(min(maxval, enr))) * .1

for cl in sorted(range(n_clusters), key=lambda v:clorder[v]):
    for i, label in enumerate(clu.labels_):
        if label != cl:
            continue
        name = motif_names[i]
        ostr = '{}\t{}\t{}'.format(name, label, motifs['Motif'][name])
        if cnv:
            cnv.setFillColor('black')
            cnv.drawString(left - 70, y, name)
            cnv.drawString(left + xstep * (len(celltypes) + .5), y, motifs['Motif'][name])
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
                xsize = enrichment2size(enrichment) * xstep#(1.0 + min(2, np.log2(enrichment))) * xstep * .15
                if enrichment <= 1.0:
                    xsize = 1.0
                    cnv.setFillColor(gray)
                    roundlevel = 0.5
                else:
                    cnv.setFillColor(colors[j])
                    roundlevel = xsize / 2
                ysize = xsize / xstep * ystep * 2
                cnv.roundRect(cx - xsize, cy - ysize, xsize * 2, ysize * 2, roundlevel, 0, fill=1)
                    #cnv.rect(cx + xstep * j - xsize, cy - ysize, xsize * 2, ysize * 2, 0, fill=1)
        y -= ystep
        if 1 or vmax > ethr: ost.write(ostr + '\n')
if args.o: ost.close()
if cnv:
    for i in range(0, 5):
        if i == 0:
            xsize = 1.0
        else:
            xsize = enrichment2size(2 ** (i - 1)) * xstep
        # logenr = i + 1
        # xsize = enrichment2size( 2 ** logenr) * xstep
        ysize = xsize / xstep * ystep * 2
        roundlevel = xsize / 2
        cnv.setFillColor('navy')
        cnv.roundRect((i  * 2)  * xstep * .5 + left - xsize, y - ysize - 100, xsize * 2, ysize * 2, roundlevel, 0, 1)
    cnv.save()