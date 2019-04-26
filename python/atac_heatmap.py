import os, sys, re, argparse
import argparse
import pandas as pd
import numpy as np
import sklearn.cluster

parser = argparse.ArgumentParser()
parser.add_argument('-i', default='THS/500bp.tsv')
parser.add_argument('-o', default='output.html')
parser.add_argument('-n', type=int, default=10, metavar='number', help='number of clusters')
parser.add_argument('--fdr', default=0.05, type=float)
parser.add_argument('--path', default='ebseq.out')
parser.add_argument('--quantile', action='store_true')
parser.add_argument('--raw', action='store_true')
parser.add_argument('--relative', action='store_true')
args = parser.parse_args()
filename_output = args.o

def quantile_normalization(matrix):
    ranks = matrix.stack().groupby(matrix.rank(method='first').stack().astype(int)).mean()
    return matrix.rank(method='min').stack().astype(int).map(ranks).unstack()


filename_input = args.i#'THS/500bp.tsv'
fdr = 0.05
displayall_mode = args.raw

paths = {}

conditions = []
levels = []
data = {}
nums = []
with open(filename_input) as fi:
    header = fi.readline().split('\t')
    for item in header:
        elems = item.split('_', 1)
        # if item == 'MUTE_rep1' or item == 'MUTE_rep2': 
        #     print('skip', item)
        #     continue
        if len(elems) >= 2:
            cnd = elems[0]
            conditions.append(cnd)
            if len(levels) == 0 or levels[-1] != cnd:
                levels.append(cnd)
                nums.append(1)
            else:
                nums[-1] += 1
    num_points = len(nums)
atac = pd.read_csv(filename_input, sep='\t', index_col=0)
genic = [x_ for x_ in atac.index if x_.find('rRNA') < 0]
if args.quantile:
    matrix = quantile_normalization(atac.loc[genic])
else:
    matrix = atac.loc[genic]
values = matrix.values
import scipy.stats

significant = []
valmat = []

for i, elem in enumerate(matrix.index):
    row = values[i]
    col = 0
    r_ = []
    valset = []
    for j in range(num_points):
        r_.append(np.mean(row[col:col+nums[j]]))
        valset.append(row[col:col+nums[j]])
        col += nums[j]
    vmin = min(r_)
    vmax = max(r_)
    pval = scipy.stats.f_oneway(*valset).pvalue
    if vmax > vmin * 2 and pval < 0.05:
        ostr = matrix.index[i]
        for v in r_:
            ostr += '\t{:.1f}'.format(v)
        # print(ostr)
        significant.append(matrix.index[i])
    else:
        g_ = matrix.index[i]
        if g_.find('(GC1)') >= 0:
            significant.append(g_)            
    data[elem] = r_

selected_matrix = pd.DataFrame([data[x_] for x_ in significant], index=significant, columns=levels)
#selected_matrix = matrix.loc[significant]#pd.DataFrame([data[x_] for x_ in significant], index=significant, columns=levels)
n_clusters = args.n
clu = sklearn.cluster.AgglomerativeClustering(n_clusters)#, affinity='cosine', linkage='complete')

relative = []
# selected_matrix = matrix
selected_matrix[selected_matrix<1] = 1
lm = np.log2(selected_matrix)
vmin, vmax = 0, 10
if args.relative:
    lm = np.subtract(lm, np.median(lm, axis=1).reshape(lm.shape[0], 1))
    vmin, vmax = -2, 2
# lm = np.subtract(lm, np.median(lm, axis=1).reshape(lm.shape[0], 1))
# lm = selected_matrix
clu.fit_predict(lm)#np.log(selected_matrix))

#clone = matrix.copy()
#clone[clone<1] = 1
if displayall_mode:
    s_ = matrix.loc[significant]
    s_[s_<1] = 1
    display_data = np.log2(s_)
    if args.relative:
        display_data = np.subtract(display_data, np.median(display_data, axis=1).reshape(display_data.shape[0], 1))
else:
    display_data = lm
# display_data = lm
# display_data = matrix.copy()
# display_data[display_data<1] = 1
# display_data =np.log2(display_data)
# display_data = np.subtract(display_data, np.median(display_data, axis=1).reshape(display_data.shape[0], 1))

def load_normalized_expression(filename):
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
            sample2column[name].append(i)
        for line in fi:
            items = line.strip().split('\t')
            gene = items[0].split('.')[0]
            if gene not in data: 
                data[gene] = np.array([float(x_) for x_ in items[1:]])
            else:
                data[gene] += np.array([float(x_) for x_ in items[1:]])
    genes = sorted(data.keys())
    df = pd.DataFrame([data[g] for g in genes], index=genes, columns=columns)
    normalized = quantile_normalization(df)
    # aggregation
    aggregated = []
    colsets = []
    for sample in samples:
        cols = []
        for col in sample2column[sample]:
            cols.append(col)
        colsets.append(cols)
    for gene in genes:
        row = normalized.loc[gene].values
        vals = []
        for cols in colsets:
            vals.append(np.mean([row[col] for col in cols]))
        aggregated.append(vals)
    return pd.DataFrame(aggregated, columns=samples, index=genes)

def load_expression(filename):
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

expression = load_normalized_expression('SRP043607.tpm')

# display_data = lm
z = []
x = display_data.columns
y = []

# for c in range(n_clusters):
#     for i, s in enumerate(significant):
#         if clu.labels_[i] == c:

expressedgenes = set(expression.index)

import scipy.stats
svalues = lm.values
clusters = [[] for i in range(n_clusters)]
n_points = len(svalues[0])
stack = np.zeros((n_clusters, n_points))
for i, l in enumerate(clu.labels_):
    clusters[l].append(i)
    stack[l] += np.array(svalues[i])
slope = np.zeros(n_clusters)
for i in range(n_clusters):
    print(stack[i])
    slope[i] = scipy.stats.linregress(np.arange(n_points), stack[i] / max(1, len(clusters[i])))[0]
print(slope)
ordered_clusters = sorted(np.arange(0, n_clusters), key=lambda x_:slope[x_], reverse=True)
print(ordered_clusters)
for cn in ordered_clusters:
    n_expr = 0    
    cluster_size =0
    expr_genes = []
    for index in clusters[cn]:
        location = display_data.index[index]
        cluster_size += 1
        ma = re.match("AT\\dG\\d+", location)
        if ma is not None:
            agi = ma.group(0)
            if agi in expressedgenes:
                expr_genes.append(agi)
                n_expr += 1

        y.append(location)
        ostr = '{}\t{}'.format(cn, location)
        display_row = display_data.loc[location].values
        z.append([max(vmin, min(vmax, x_)) for x_ in display_row])
        for v in display_row:
            ostr += '\t{:.1f}'.format(v)
    if n_expr > 0:
        emat = expression.loc[expr_genes]
        print(cluster_size, n_expr, [int(x*100) for x in np.mean(emat, axis=0)])
import plotly.graph_objs as go
import plotly, plotly.offline
fig = go.Figure([go.Heatmap(x=x, y=y, z=z)])
plotly.offline.plot(fig, filename=filename_output)
