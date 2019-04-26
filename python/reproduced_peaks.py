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

#if __name__ == '__main__':

if __name__ == '__main__':
    srcdir = dstdir = 'THS'
    samples = {}
    for fn in os.listdir(srcdir):
        # ('MUTE_rep4_summits.bed', None)
        # ('MUTE_rep3_peaks.xls', None)
        m = re.match('(.+?)_rep\\d+_peaks\\.xls$', fn)
        if m:
            name = m.group(1)
            if name not in samples: samples[name] = []
            samples[name].append(os.path.join(srcdir, fn))
    for name, fns in samples.items():
        with open(os.path.join(dstdir, '{}.bed'.format(name)), 'w') as fo:
            detectedpeaks = []
            chromosomes = set()
            for fn in fns:
                p = load_bed(fn)
                detectedpeaks.append(p)
                for c in p.keys(): chromosomes.add(c)
            num = 1
            for c in chromosomes:
                chpeaks = [bed.get(c, ()) for bed in detectedpeaks]
                reproduced = detect_overlap(chpeaks, 3)
                for start, stop in reproduced:
                    fo.write('{}\t{}\t{}\t{}_{}\n'.format(c, start, stop, name, num))
                    num += 1
        sys.stderr.write('{}\t{}\n'.format(name, num))

def test2():
    beddir = 'macs'
    dstdir = 'THS'
    with open('description.txt') as fi:
        for line in fi:
            items = line.strip().split('\t', 1)
            if len(items) != 2: continue
            name = items[0]
            samples = [x_.strip() for x_ in items[1].split(',')]
            peakfiles = []
            for sample in samples:
                bed = os.path.join(beddir, sample + '_peaks.xls')
                if os.path.exists(bed) is False:
                    sys.stderr.write('{}:{} not found\n'.format(name, bed))
                else:
                    peakfiles.append(load_bed(bed))
            if len(peakfiles) == 0:
                sys.stderr.write('{} has no peaks\n'.format(name))
                continue
            n = 0
            with open(os.path.join(dstdir, name + '.bed'), 'w') as fo:
                fo.write('track type=bed name="{} reproduced THS"\n'.format(name))
                for chrom in peakfiles[0].keys():
                    peaks = [bed.get(chrom, ()) for bed in peakfiles]
                    reproduced = detect_overlap(peaks, 2)
                    for start, stop in reproduced:
                        fo.write('{}\t{}\t{}\t{}_{}\n'.format(chrom, start, stop, name, n))
                        n += 1
            sys.stderr.write('{} : {} peaks output\n'.format(name, n))
                    

def test():
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
