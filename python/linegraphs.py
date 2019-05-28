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
args = parser.parse_args()
vmax = args.vmax

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


