import os, sys, subprocess
import tempfile

filename = 'atac_500bp.tsv'
filename = sys.argv[1]
states = {}
filename_output = filename + '.out'
filename_annotation = 'annotation.tsv'
filename_script = 'script.R'

with open(filename) as fi, open(filename_annotation, 'w') as fo:
    fo.write('\tSample\tCondition\tTime\tBatch\n')
    header = fi.readline().split('\t')
    # print(header)
    timepoints = []
    states = {}
    for item in header[1:]:
        item = item.strip()
        state = item.split('_')[0]
        if state not in states:
            states[state] = len(states) + 1
        tm = states[state]
        fo.write('{0}\t{0}\tcase\t{1}\tNULL\n'.format(item, tm))

"""
Sample Condition Time Batch
"""

script = """
suppressPackageStartupMessages(library(ImpulseDE2))
suppressPackageStartupMessages(library(splines))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(BiocParallel))

## load
countdata <- read.csv("{}", sep="\\t", header=T, row.names=1)
annotation <- read.csv("{}", sep="\\t", header=T, row.names=1)
## simulate data
lsSim <- simulateDataSetImpulseDE2(
    vecTimePointsA=rep(1:6, 3),
    vecTimePointsB=NULL,
    vecBatchesA=rep("B_NULL", 6*3),
    vecBatchesB=NULL,
    scaNConst=100,
    scaNImp=50,
    scaNLin=50,
    scaNSig=50,
    scaNRand=0,
    scaMumax=1000,
    scaSeedInit = 1)
## run DESeq2
print(annotation)
splines <- ns(annotation$Time, df=4)
colnames(splines) <- paste0("spline", seq(1, dim(splines)[2]))
print(splines)

# create DESeq2 object with spline basis vectors as individual linear predictors
dds <- DESeqDataSetFromMatrix(
    countData = countdata,
    colData = splines,
    design = ~spline1 + spline2 + spline3 + spline4)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dds <- estimateDispersionsFit(dds)

dds <- estimateDispersionsMAP(dds, outlierSD = 10)
# perform a log-likelihood ratio test
dds <- nbinomLRT(dds, full = ~spline1 + spline2 + spline3 + spline4,
                 reduced = ~1)
write.table(results(dds), "{}", sep="\t", col.names=T, quote=F)
""".format(filename, filename_annotation, filename_output)

with open(filename_script, 'w') as fo:
    fo.write(script)
subprocess.Popen(['Rscript', filename_script]).wait()

os.unlink(filename_script)
os.unlink(filename_annotation)