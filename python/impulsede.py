import os, sys, subprocess, re
import tempfile
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', default=None)
parser.add_argument('-o', default=None)
# parser.add_argument('-a', default='annotation.tsv')
args = parser.parse_args()

# filename = 'atac_500bp.tsv'
filename = sys.argv[1]
states = {}
filename_output = filename + '.out'
filename_annotation = 'annotation.tsv'
filename_script = 'script.R'

filename = args.i
filename_output = args.o if args.o is not None else filename + '.out'
filename_annotation = tempfile.mktemp('.tsv')
filename_script = tempfile.mktemp('.R')

with open(filename) as fi, open(filename_annotation, 'w') as fo:
    fo.write('ID\tSample\tCondition\tTime\tBatch\n')
    header = fi.readline().split('\t')
    # print(header)
    timepoints = []
    states = {}
    for item in header[1:]:
        item = item.strip()
        state, rep = item.split('_', 1)
        rep = re.sub('rep', '', rep)
        if state not in states:
            states[state] = len(states) + 1
        tm = states[state]
        fo.write('{0}\t{0}\tcase\t{1}\tB_{2}\n'.format(item, tm, rep))
        #fo.write('{0}\t{0}\tcase\t{1}\tB_NULL\n'.format(item, tm, rep))

"""
Sample Condition Time Batch
"""

script = """
suppressPackageStartupMessages(library(ImpulseDE2))
suppressPackageStartupMessages(library(splines))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(BiocParallel))

## load
countdata <- read.csv("{0}", sep="\\t", header=T, row.names=1)
annotation <- read.csv("{1}", sep="\\t", header=T, row.names=1)
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
write.table(results(dds), "{2}", sep="\t", col.names=T, quote=F)
""".format(filename, filename_annotation, filename_output)


simulation = """
lsSimulatedData <- simulateDataSetImpulseDE2(
  vecTimePointsA   = rep(seq(1,8),3),
  vecTimePointsB   = NULL,
  vecBatchesA      = c(rep("B1",8), rep("B2",8), rep("B3",8)),
  vecBatchesB      = NULL,
  scaNConst        = 0,
  scaNImp          = 100,
  scaNLin          = 0,
  scaNSig          = 0,
  scaMuBatchEffect = 1,
  scaSDBatchEffect = 0.2,
  dirOutSimulation = NULL)

countdata <- lsSimulatedData$matObservedCounts
annotation <- lsSimulatedData$dfAnnotation
head(countdata)
dim(countdata)
dim(annotation)
"""

script = """
suppressPackageStartupMessages(library(ImpulseDE2))
suppressPackageStartupMessages(library(splines))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(BiocParallel))





countdata <- read.csv("{0}", sep="\\t", header=T, row.names=1)
annotation <- read.csv("{1}", sep="\\t", header=T, row.names=1)


for (col in colnames(countdata)) {{
  countdata[col] <- as.numeric(unlist(countdata[col]))
}}
countdata <- as.matrix(countdata)

annotation$Time<-as.numeric(unlist(annotation$Time))

sapply(countdata, class)
sapply(annotation, class)

head(countdata)

print(annotation)

## lsSimulatedData$matObservedCounts, 
## lsSimulatedData$dfAnnotation
objectImpulseDE2 <- runImpulseDE2(
  matCountData           = countdata,
  dfAnnotation           = annotation,
  boolCaseCtrl           = FALSE,
  vecConfounders         = c("Batch"),
  boolIdentifyTransients = TRUE,
  scaNProc               = 1 )

head(objectImpulseDE2$dfImpulseDE2Results)
write.table(objectImpulseDE2$dfImpulseDE2Results, "{2}", sep="\t", col.names=T, quote=F)

""".format(filename, filename_annotation, filename_output)

with open(filename_script, 'w') as fo:
    # print(script)
    fo.write(script)
# print('#CHECK')
# with open(filename) as fi:
#     n = len(fi.readline().strip().split('\t'))
#     for line in fi:
#         items = line.strip().split('\t')
#         if len(items) != n:
#             sys.stderr.write(line)
#         for i, v in enumerate(items[1:]):
#             try:
#                 v = int(v)
#             except:
#                 sys.stderr.write(line + '\n')
#                 break
# exit()

subprocess.Popen(['Rscript', filename_script]).wait()

os.unlink(filename_script)
os.unlink(filename_annotation)