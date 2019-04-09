import os, sys, re

def detect_overlap(peakdata, threshold=0):
    if threshold == 0:
        threshold = len(peakdata)
    peaks = []
    for data in peakdata:
        prev = -1
        for start, stop in data:
            if start < prev: start = prev
            peaks.append((start, stop))
    peaks = sorted(peaks, key=lambda x:x[0])
    i = 0
    overlapped = []
    while i < len(peaks):
        pi = peaks[i]
        j = i + 1
        num = 1
        peak_start, peak_stop = pi
        while j < len(peaks):
            pj = peaks[j]
            if peak_stop < pj[0]:
                break
            peak_start = max(peak_start, pj[0])
            peak_stop = min(peak_stop, pj[1])
            num += 1
            if num >= threshold:
                # print('accept {}-{}:{}, {}'.format(peak_start, peak_stop, num, str(peaks[i:j+1])))
                overlapped.append((peak_start, peak_stop))
            j += 1
        i += 1
    merged = []
    i = 0
    while i < len(overlapped):
        p0 = overlapped[i]
        start = p0[0]
        stop = p0[1]
        j = i + 1
        while j < len(overlapped):
            pj = overlapped[j]
            if pj[0] <= stop:
                stop = pj[1]
            else:
                break
            j += 1
        merged.append((start, stop))
        i += 1
    return merged

def load_bed(filename):
    data = {}
    with open(filename) as fi:
        for line in fi:
            items = line.strip().split('\t')
            if len(items) >= 3 and items[1].isdigit():
                chrom = items[0]
                start = int(items[1])
                stop = int(items[2]) + 1
                if chrom not in data:
                    data[chrom] = []
                data[chrom].append((start, stop))
    return data

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', nargs='+')
parser.add_argument('-o', default='out.bed')
#parser.add_argument('-s', nargs='+', default=None)
args = parser.parse_args()

bedfiles = []
for fn in args.i:
    bedfiles.append(load_bed(fn))
results = {}
with open(args.o, 'w') as fo:
    name = os.path.basename(args.o).split('.')[0]
    n = 1
    for chrom in bedfiles[0].keys():
        peaks = [bed.get(chrom, ()) for bed in bedfiles]
        reproduced = detect_overlap(peaks, 2)
        for start, stop in reproduced:
            fo.write('{}\t{}\t{}\t{}_{}\n'.format(chrom, start, stop, name, n))
            n += 1
        results[chrom] = reproduced
# if args.s is not None:
#     summits = [load_bed(fn) for fn in args.s]
#     for chrom in results.keys():
#         peaks = results[chrom]



    
