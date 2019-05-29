"""
Read results.tsv file which is generated by atac_heatmap.py and draw graphs

"""

import os, sys, re
import reportlab.pdfgen.canvas
import tkgraph
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i', metavar='filename', default='MyData/heatmap/results.tsv')
parser.add_argument('-o', default='out.pdf')
parser.add_argument('--vmax', default=3.0, type=float)
parser.add_argument('--expression', action='store_true')
parser.add_argument('--summary', default=None)
args = parser.parse_args()
vmax = args.vmax
filename_summary = args.summary
if filename_summary is not None:
    ordered_clusters = []
    with open(filename_summary) as fi:
        fi.readline()
        for line in fi:
            items = line.strip().split('\t')    
            if items[0].isdigit():
                ordered_clusters.append(int(items[0]))
    ordered_clusters = list(reversed(ordered_clusters))
else:
    ordered_clusters = None
data = []
clusters = []
#locations = 
with open(args.i) as fi:
    fi.readline()
    for line in fi:
        items = line.strip().split('\t')
        cluster = int(items[0])
        if cluster >= len(data):
            # print(len(data), cluster)
            for i in range(len(data), cluster + 1):
                data.append([])
                clusters.append([])
        gene = items[1]
        location = items[2]
        clusters[cluster].append(gene)
        row = []
        for item in items[3:]:
            if item == '.' or item == 'na':
                row.append(None)
            else:
                row.append(float(item))
        data[cluster].append(row)
        # print(row)
num_clusters = len(data)
import scipy.stats

atac_columns = [[3,4,5], [6,7,8], [9,10,11,12,13], [14,15,16], [17,18,19]]
expr_columns = [20, 21, 22, 23, 24]
if ordered_clusters is None:
    ordered_clusters = list(range(num_clusters))

import reportlab.pdfgen.canvas
import tkgraph
gr = tkgraph.ScatterDiagram('scatters.pdf', tkgraph.Axis(1,3,50,150), tkgraph.Axis(-1,3,600,700))
gr.set_color(0,0, 128)
gr.set_dot_size(3)

# cnv = reportlab.pdfgen.canvas.Canvas('correlation_graph.pdf')
row = 0
for cluster in ordered_clusters:
    data_of_cluster = data[cluster]
    edata = [[] for i in range(5)]
    adata = [[] for i in range(5)]
    nums = [0] * 5
    for vals in data_of_cluster:
        n = 0
        for i in range(5):
            xs = [vals[col - 3] for col in atac_columns[i]]
            x = np.mean(xs)
            y = vals[-5 + i]
            if x is not None and y is not None:
                n += 1
                adata[i].append(x + 1)
                edata[i].append(y + .01)
    # for i in range(5):
    #     if nums[i] > 0:
    #         adata[i] /= nums[i]
    #         edata[i] /= nums[i]
        # for cols in rangeatac_columns:
        #     # try:
        #     atac.append(np.mean([(datum[x_ - 3]) for x_ in cols]))
        #     # except:
        #     #     # print(datum, cols)
        #     #     raise
        # expr = []
        # for col in expr_columns:
        #     val = datum[col - 3]
        #     if val is None:
        #         expr = None
        #         break
        #     else:
        #         expr.append(val + .1)
        # if expr is not None:
        #     adata.append(atac)
        #     edata.append(expr)
    # print(len(adata), len(edata))
    # if len(adata) > 1 and len(edata) > 1:
    ostr = '{}\t{}'.format(cluster, len(adata[0]))
    available = False
    for i in range(5):
        if len(adata[i]) < 3:
            ostr += '\t.; .; .'
        else:
            avals = adata[i]
            evals = edata[i]
            if len(avals) > 0:
                if i == 0:
                    gr.drawString(gr.left, gr.top + 5, str(i))
                for j in range(len(avals)):
                    x = np.log10(avals[j] + 1)
                    y = np.log10(evals[j] + .1)
                    gr.set_point(x, y)
                gr.shift(xspacer=10, yspacer=10)
            r, p = scipy.stats.pearsonr(np.log(avals), np.log(evals))#adata[i], edata[i])
            # ostr += '\t{:.2f}; {:.2f}; {:.3f}'.format(np.mean(avals), np.mean(evals), r)
            ostr += '\t{:.4f}'.format(r)#np.mean(avals), np.mean(evals), r)
            available = True
    if available:
        # for i, r in enumerate(rvals):

        print(ostr)
        row += 1
    #     avalues = np.array([adatum[i] for adatum in adata])
    #     evalues = np.array([edatum[i] for edatum in edata])
    #     amean = np.mean(avalues)
    #     evalues[evalues<0.1] = 0.1
    #     emean = np.mean(np.log(evalues))
    #     ostr += '\t{:.3f}; {:.3f}; {:.3f}'.format(amean, emean, r)
    # print(ostr)
gr.save()
exit()






            
        


data = []
with open(args.i) as fi:
    fi.readline()
    for line in fi:
        items = line.strip().split('\t')
        cluster = int(items[0])
        if cluster >= len(data):
            # print(len(data), cluster)
            for i in range(len(data), cluster + 1):
                data.append([])
        gene = items[1]
        location = items[2]
        if args.expression:
            elems = items[-5:]
            if '.' in elems or 'nan' in elems:
                continue
            values = np.array([float(x_) for x_ in elems])
            # values = np.array([max(0.01, float(x_)) for x_ in elems])
            #values = np.log(values) 

            mean = np.mean(values)
            if mean < 0.01: continue
            # if mean < 1.0: continue
            # print(mean, values)

            sd = np.std(values)
            #values = np.log2(values) - np.log2(np.median(values))
            values = (values - mean) / sd
            # print(values)
        else:
            values = []
            for cols in [[3,4,5], [6,7,8], [9,10,11,12,13], [14,15,16], [17,18,19]]:
                s = 0
                n = len(cols)
                for col in cols:
                    s += float(items[col])
                values.append(s / n)
        data[cluster].append(values)

gr = tkgraph.Graph(args.o, tkgraph.Axis(0, 5, 50, 175), tkgraph.Axis(-vmax, vmax, 600, 750, major_tick=vmax))
gr.setFont('Helvetica', 8)
def generate_path(graph=gr, values=values, path=None):
    if path is None:
        path = gr.beginPath()
    for i, v in enumerate(values):
        if np.isnan(v) or np.isinf(v): 
            v = 0
    # if not np.isinf(v) and not np.isnan(v):
        x = gr.x(i + .5)
        y = gr.y(v)
        # print(v, x, y)
        if i == 0:
            path.moveTo(x, y)
        else:
            path.lineTo(x, y)
    return path
for i, valset in enumerate(data):
    path = None
    accum = [[] for i in range(5)]#[0] * 5
    for values in valset:
        for j, v in enumerate(values):
            accum[j].append(v)# += v
        path = generate_path(graph=gr, values=values, path=path)
    # print(i, len(valset))
    if path is not None:
        gr.setStrokeColorRGB(.8, .8, .8)
        gr.setLineWidth(.1)
        gr.drawPath(path)
        # path = generate_path(graph=gr, values=[x_ / len(valset) for x_ in accum])
        repr = [np.mean(x_) for x_ in accum]
        print(i, repr)
        path = generate_path(graph=gr, values=repr)#np.median(accum, axis=0))#[np.medianx_ / len(valset) for x_ in accum])
        gr.setStrokeColor('crimson')
        gr.setLineWidth(1)
        gr.drawPath(path)
        gr.drawString(gr.right - 50, gr.top - 20, '{} genes'.format(len(valset)))
        gr.shift(yspacer=10, xspacer=20)
gr.save()

