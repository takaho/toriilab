import os, sys, re, argparse
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i', default='THS/500bp.tsv')
parser.add_argument('--fdr', default=0.05, type=float)
parser.add_argument('--path', default='ebseq.out')
parser.add_argument('--quantile', action='store_true')
args = parser.parse_args()

def quantile_normalization(matrix):
    ranks = matrix.stack().groupby(matrix.rank(method='first').stack().astype(int)).mean()
    return matrix.rank(method='min').stack().astype(int).map(ranks).unstack()

# import pandas as pd
# mat = pd.read_csv('dummy.tsv', sep='\t', index_col=0)
# print(mat)
# print(quantile_normalization(mat))
# exit()

filename_input = args.i#'THS/500bp.tsv'
filename_output = args.path#'test.out'
fdr = 0.05


paths = {}

conditions = []
levels = []
data = {}
nums = []
with open(filename_input) as fi:
    header = fi.readline().split('\t')
    for item in header:
        elems = item.split('_', 1)
        if len(elems) >= 2:
            cnd = elems[0]
            conditions.append(cnd)
            if len(levels) == 0 or levels[-1] != cnd:
                levels.append(cnd)
                nums.append(1)
            else:
                nums[-1] += 1
    num_points = len(nums)
matrix = quantile_normalization(pd.read_csv(filename_input, sep='\t', index_col=0))
values = matrix.values
import scipy.stats

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
    if vmax > vmin * 2 and pval < 0.01:
        ostr = matrix.index[i]
        for v in r_:
            ostr += '\t{:.1f}'.format(v)
        print(ostr)
    data[elem] = r_

    # for line in fi:
    #     items = line.strip().split('\t')
    #     values = []#float(x) for x in items[1:]]
    #     col = 1
    #     for i in range(num_points):
    #         n_ = nums[i]
    #         v_ = 0
    #         for j in range(n_):
    #             v_ += float(items[col])
    #             col += 1
    #         values.append(v_ / n_)
    #     data[items[0]] = values
exit()

if os.path.exists(filename_output) is False:

    def convert2Rarray(values):
        out = ''
        for item in values:
            if out != '': out += ','
            out += '"{}"'.format(item)
        return out

    script = """
library(EBSeqHMM)
data(GeneExampleData)
df <- read.csv("{inputfilename}", sep="\\t", row.names=1)
countdata <- as.matrix(df)

Conditions <- factor(c({conditions}), levels=c({levels}))
Sizes <- MedianNorm(countdata)
EBSeqHMMGeneOut <- EBSeqHMMTest(Data=countdata, sizeFactors=Sizes, Conditions=Conditions, UpdateRd=2)
GeneDECalls <- GetDECalls(EBSeqHMMGeneOut, FDR={fdr})
write.table()
write.table(GeneDECallas, filename="{outputfilename}", sep="\\t", quote=FALSE)
    """.format(inputfilename=filename_input, conditions=convert2Rarray(conditions), levels=convert2Rarray(levels), fdr=fdr, outputfilename=filename_output)

    import tempfile

    ft = tempfile.mktemp('.R')
    with open(ft, 'w') as fo:
        fo.write(script)

    cmd = 'Rscript', script
    subprocess.Popen(cmd).wait()
    with open(filename_output) as fi:
        fi.readline()
        for line in fi:
            items = line.strip().split('\t')
            paths[items[0]] = items[1], float(items[2])
    os.unlink(ft)
else:
    with open(filename_output) as fi:
        fi.readline()
        for line in fi:
            items = line.strip().split('\t')
            paths[items[0]] = items[1], float(items[2])


patterns = set([x_[0] for x_ in paths.values()])
locations = sorted(data.keys())
for pattern in sorted(patterns):
    for loc in locations:
        if loc in paths and paths[loc][0] == pattern:
            values = data[loc]
            vmax = max(values)
            vmin = min(values)
            if vmax > vmin * 2:
                ostr = '{}\t{}'.format(pattern, loc)
                for v in values:
                    ostr += '\t{:.2f}'.format(v)
                print(ostr)