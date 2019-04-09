import os, sys, re

def load_peaks(filename):
    genes = set()
    with open(filename) as fi:
        for line in fi:
            items = line.strip().split('\t')
            for elem in items[4].split('//'):
                if elem.startswith('TSS_'):
                    gene = elem.split('_', 1)[1].split('.')[0]
                    genes.add(gene)
    return genes
def load_tpm(filename):
    columns = [[1,2], [3, 4, 5, ], [6,7], [8,9], [10,11]]
    #columns = [[1,2], [3,4,5], [6,7], [8,9], [10,11]]
    matrix = {}
    with open(filename) as fi:
        names = fi.readline().strip().split('\t')[1:]
        genes = []
        for line in fi:
            items = line.strip().split('\t')
            gene = items[0]
            genes.append(gene)
            values = []
            for cols in columns:
                a = 0
                for col in cols:
                    a += float(items[col])
                values.append(a / len(cols))
            matrix[gene] = values
            #matrix.append(values)
    return matrix

import numpy as np                         
labels = 'AtML1', 'SPCH', 'MUTE', 'FAMA', 'GC1'
geneset = []
gene2group = {}
expr = load_tpm('SRP032607.tpm')
notexpressed = set()
for g, row in expr.items():
    if np.sum(row) <= 0.0:
        notexpressed.add(g)
for i, label in enumerate(labels):
    fn = 'peakcall/{}.assigned.txt'.format(label)
    genes = load_peaks(fn)
    geneset.append(genes)
    plus, minus = [], []
    for g, row in expr.items():
        #if g in notexpressed: continue
        if g in genes:
            plus.append(row[i])
        else:
            minus.append(row[i])
    lp = [np.log2(x_) for x_ in plus if x_ > 0]
    lm = [np.log2(x_) for x_ in minus if x_ > 0]
    print('{}\t{}\t{}\t{:.3f}\t{:.3f}\t{:.1f}\t{:.1f}'.format(label, len(minus), len(plus), np.median(minus), np.median(plus), np.mean(lm), np.mean(lp)))
    flag = 1 << i
    for gene in genes:
        gene2group[gene] = gene2group.get(gene, 0) + flag
    

accum = [0] * len(labels)
for row in expr.values():
    for i, v in enumerate(row):
        accum[i] += v
averages = [a / len(expr) for a in accum]

groups = [[] for i in range(max(gene2group.values()) + 1)]
for gene, gr in gene2group.items():
    groups[gr].append(gene)
num_steps = len(labels)


for i in range(max(gene2group.values()) + 1):
    n = len(groups[i])
    if n < 100: continue
    #accum = [0] * num_steps
    subdata = [[] for j in range(num_steps)]
    vstr = ''
    for gene in groups[i]:
        if gene not in expr: 
            sys.stderr.write('{}\n'.format(gene))
            continue
        row = expr[gene]
        #print(row)
        #for index in range(len(labels)):
        for j, v in enumerate(row):
            subdata[j].append(v)
            #accum[j] += v
    #print(accum)
    #avg = [x_ / n for x_ in accum]
    for j in range(num_steps):
        #for j, a in enumerate(accum):
        if j > 0: vstr += ', '
        vstr += '{:.2f}'.format(np.mean(subdata[j]))
        #vstr += '{:.2f}'.format(np.median(subdata[j]))
    flag = ''
    for j in range(len(labels)):
        if (i & (1 << j)) != 0:
            flag += 'O'
        else:
            flag += '.'
    print('{}\t{}\t{}\t{}'.format(i, flag, len(groups[i]), vstr))
