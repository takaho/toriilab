import os, sys, re, argparse

def load_expression(filename):
    # import pandas as pd
    # import collections
    # table = pd.read_csv(filename, sep='\t', index_col=0)
    # columns = table.columns
    # samples = collections.OrderedDict()
    # for i, col in enumerate(columns):
    #     name = col.split('_')[0]
    #     if name not in samples:
    #         samples[name] = []
    #     samples[name].append(col)
    # aggregated = pd.join([np.mean(table[cols]) for cols in samples.values()])
    # aggregated.columns = samples.keys()
    # return aggregated

    data = {}
    with open(filename) as fi:
        columns = fi.readline().strip().split('\t')[1:]
        samples = []
        sample2column = {}
        for i, col in enumerate(columns):
            name = col.split('_', 1)[0]
            if len(samples) == 0 or samples[-1] != name:
                samples.append(name)
                sample2column[name] = []
            sample2column[name].append(i + 1)
        for line in fi:
            items = line.strip().split('\t')
            gene = items[0].split('.')[0]
            if gene not in data: data[gene] = [0] * len(samples)
            row = data[gene]
            for i, name in enumerate(samples):
                cols = sample2column[name]
                val = 0
                for col in cols:
                    val += float(items[col])
                row[i] += val / len(cols)
    return data
            

parser = argparse.ArgumentParser()
parser.add_argument('-i', default='THS/500bp.imp', help='InpulseDE2 output')
parser.add_argument('-t', default='THS/500bp.tpm', help='TPM')
parser.add_argument('-o', default='out.html', help='Heatmap file')
parser.add_argument('-n', type=int, default=10, metavar='number', help='number of clusters')
parser.add_argument('--colorscale', default='0,50', help='heatmap min/max')
parser.add_argument('--log', action='store_true')
parser.add_argument('--verbose', action='store_true')
parser.add_argument('--relative', action='store_true')
parser.add_argument('--qvalue', default=0.05, type=float, help='qValue threshold')
# parser.add_argument('--logfc', default=0.0, type=float, help='log2 fold change threshold')
parser.add_argument('--min-tpm', default=0.0, type=float, help='mean TPM threshold')
parser.add_argument('--known', action='store_true')
parser.add_argument('--display-all', action='store_true')
args = parser.parse_args()

verbose = args.verbose
use_log = args.log
filename_output = args.o
relative_mode = args.relative
bundle_mode = not args.display_all
use_known = args.known
# filename_count = 'THS/500bp_wt.tsv'
# filename_tpm = 'THS/500bp_wt.tpm'
# filename_imp = 'THS/500bp_wt.tsv.out'

filename_tpm = args.t
filename_imp = args.i
n_clusters = args.n

vmin, vmax = [float(x_) for x_ in args.colorscale.split(',')[0:2]]
# check

if 0:
    import scipy.cluster.hierarchy
    import pandas as pd
    filename_tpm = 'THS/500bp.tpm'
    t = pd.read_csv(filename_tpm, sep='\t', index_col=0)
    linkage = scipy.cluster.hierarchy.linkage(t.T)
    print(linkage)
    exit()

# vmax = 0.5

# filename_count = 'THS/500bp_wt.tsv'
# filename_tpm = 'THS/500bp.tpm'
# filename_imp = 'THS/500bp.tsv.out'

# filename_count = 'THS/500bp_wt.tsv'
# filename_tpm = 'THS/1k.tpm'
# filename_imp = 'THS/1k.tsv.out'
#vmax = 50

#qthr = 0.05
qthr = args.qvalue
# lthr = args.logfc
#vmax = 7
#vmin = 0.0
locations = {}
tpm_threshold = args.min_tpm
with open(filename_imp) as fi:
    fi.readline()
    for line in fi:
        if 'NA' in line: continue
        items = line.strip().split('\t')
        # print(items)
        location = items[0]
        pval = float(items[2])
        qval = float(items[3])
        mean = float(items[8])
        # values = [float(x) for x in items[1:]]
        # mean, lfc, lfcse, stat, pval, qval = values
        values = [mean, pval, qval]
        if mean > tpm_threshold and qval < qthr and (not use_known or location.find('(') > 0):
            locations[location] = values

#Kruskal
if 0:
    import scipy.stats
    N = 0
    n2 = 0
    significant = []
    with open(filename_tpm) as fi:
        columns = fi.readline().strip().split('\t')[1:]
        cellcolumns = []
        celltypes = []
        for i, c in enumerate(columns):
            celltype = c.split('_', 1)[0]
            integrated = False
            for j, ct in enumerate(celltypes):
                if ct == celltype:
                    cellcolumns[j].append(i)
                    integrated = True
                    break
            if not integrated:
                celltypes.append(celltype)
                cellcolumns.append([])
                cellcolumns[-1].append(i)
        for line in fi:
            items = line.strip().split('\t')
            gene = items[0]
            valset = []
            for cols in cellcolumns:
                valset.append([float(items[c+1]) for c in cols])
            try:
                # pval = scipy.stats.kruskal(*valset).pvalue
                pval = scipy.stats.f_oneway(*valset).pvalue
            except:
                pval = 1.0
            if pval <0.01:
                # print(gene, valset)
                if gene in locations:
                    n2 += 1
                significant.append(gene)
            N += 1
    print(significant)
    print(locations.keys())
    n0 = len(locations)
    n1 = len(significant)
    #n2 = len([g for g in significant if g in locations])
    print(n0, n1, n2, N - n0 - n1 + n2)
    locations = set(significant)
genes = []
matrix = []

import numpy as np
celltypes = []
cellcolumns = []
cellvalues = []
# vmin = -vmax
# vmin, vmax = 0, 5
# vmin, vmax = -1, 2
with open(filename_tpm) as fi:
    columns = fi.readline().strip().split('\t')[1:]
    cellcolumns = []
    for i, c in enumerate(columns):
        celltype = c.split('_', 1)[0]
        integrated = False
        for j, ct in enumerate(celltypes):
            if ct == celltype:
                cellcolumns[j].append(i)
                integrated = True
                break
        if not integrated:
            celltypes.append(celltype)
            cellcolumns.append([])
            cellcolumns[-1].append(i)
    print(celltypes, cellcolumns)
    for line in fi:
        items = line.split('\t')
        if items[0] not in locations: continue
        # mean, lfc, lfcse, stat, pval, qval = values = locations[items[0]]
        gene, position = items[0].split(';', 1)
        if gene not in genes:
            values = [float(x) for x in items[1:]]
            med = np.median(values)
            if use_log:
                logvalues = [np.log2(max(.0001, x)) for x in values]
                values = logvalues
                if relative_mode:
                    values = [(x_ - np.log2(med)) for x_ in values]
            elif relative_mode:
                values = [np.log2(x_ / med) for x_ in values]
            
            # logratio = [np.log2(max(.1, x) / med) for x in values]
            # logvalues = [np.log2(max(.0001, x)) for x in values]
            # # values = [np.log2(max(.001, v)) for v in values]
            # values = logvalues
            # matrix.append(logratio)
            matrix.append(values)#[np.log2(max(.001, v)) for v in values])
            aggregated = []
            for cols in cellcolumns:
                s = 0
                for c in cols:
                    s += values[c]
                aggregated.append(s / len(cols))
            # cellvalues.append([np.log2(v/med) for v in aggregated])
            cellvalues.append(aggregated)
            genes.append(gene)

import sklearn.cluster
import pandas as pd

df = pd.DataFrame(matrix, index=genes, columns=columns)
avg = pd.DataFrame(cellvalues, index=genes, columns=celltypes)

if bundle_mode:
    heatmapdata = avg
else:
    heatmapdata = df

clu = sklearn.cluster.AgglomerativeClustering(n_clusters=n_clusters)
clu.fit_predict(avg)
import plotly
import plotly.graph_objs as go

z = []
profiles = []
# vmin, vmax = 1, 6

matrix = heatmapdata.values
nums = [0] * n_clusters
score = [0] * n_clusters
for cn in range(n_clusters):
    for i in range(len(genes)):
        l = clu.labels_[i]
        score[l] += matrix[i][0] - matrix[i][-1]
        nums[l] += 1
ordered_clusters = list(sorted(range(n_clusters), key=lambda cn:score[cn] / max(1, nums[cn])))

expression = load_expression('SRP043607.tpm')

ordered_genes = []
for cn in ordered_clusters:
    accum = [0] * df.shape[1]
    n = 0
    trn = None
    num_expr = 0
    for i, g in enumerate(genes):
        if clu.labels_[i] == cn:
            m = re.match('AT\\dG\\d+', g)
            if m:
                expr = expression.get(m.group(0), (0, ))
                if sum(expr) > 1:
                    if trn is None:
                        trn = np.array(expr)
                    else:
                        trn += np.array(expr)
                    num_expr += 1
            ost = '{}\t{}'.format(cn, g)
            z.append([max(vmin,min(vmax, x_)) for x_ in heatmapdata.loc[g].values])
            for j, v in enumerate(heatmapdata.loc[g].values):
                accum[j] += v
                ost += '\t{:.2f}'.format(v)
            ordered_genes.append(g)
            print(ost)
            n += 1
    profiles = [x_ / max(1, n) for x_ in accum]
    if num_expr > 0:
        trn /= num_expr
        print(cn, str(trn))
print(profiles)
hm = go.Heatmap(z=z, x=heatmapdata.columns, y=ordered_genes)
fig = go.Figure([hm, ])
plotly.offline.plot(fig, filename=filename_output)
