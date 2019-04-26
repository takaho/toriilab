import os, sys, re, subprocess
import tempfile
import datetime

l2s = {}
s2l = {}
with open('/mnt/smb/tae/torii/nextseq/description.txt') as fi:
    for line in fi:
            #items = re.split('\\s+', line.strip(), 1)
        items = line.strip().split('\t')
        # print(items)
        if len(items) >= 2:
            name = items[0]
            labels = [x.strip() for x in items[1].split(',')]
            if name not in s2l:
                s2l[name] = []
            slot = len(s2l[name])
            for label in labels:
                l2s[label] = name, slot
            s2l[name].append(labels)
#            s2l[name] = labels
        
bamdir = 'bamfiles'
#bamdir  = '.'
bamfiles = {}
for sample in s2l.keys(): 
    slots = s2l[sample]
    for i in range(len(slots)):
        name = '{}_rep{}'.format(sample, i + 1)
        fns = []
        for label in slots[i]:
            fn = os.path.join(bamdir, label + '.bam')
            if os.path.exists(fn) and os.path.getsize(fn) > 1000000:
                fns.append(fn)
            else:
                sys.stderr.write('{} was not found\n'.format(fn))
        bamfiles[name] = fns

genomesize = {}
centromere = []
dstdir = 'macs2'
dstdir = 'peakcall'
dstdir = 'THS'

if not os.path.exists(dstdir): os.makedirs(dstdir)

nuc_len = 0
for name, fns in bamfiles.items():
    #print(name)
    fn_out = os.path.join(dstdir, name + '_peaks.xls')
    if os.path.exists(fn_out):# and os.path.getsize(fn_out) > 100:
        sys.stderr.write('skip {} \n'.format(name))
        continue
    sys.stderr.write('[{}] {}\n'.format(datetime.datetime.now().strftime('%H:%M:%S'), name))
    if len(fns) == 0: continue
    if len(genomesize) == 0:
        # cmd = 'samtools', 'view', '-H', fns[0]
        # proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        # for line in proc.stdout:
        #     m = re.match('@SQ\\s+SN:(\\S+)\\s+LN:(\\d+)', line.decode('ascii'))
        #     if m:
        #         genomesize[m.group(1)] = int(m.group(2))
        #         if re.match('chr\\d$', m.group(1)):
        #             nuc_len += int(m.group(2))
        # proc.stdout.close()
        # proc.wait()
        with open('TAIR10_nocentromere.bed') as fi:
            gs_ = {}
            for line in fi:
                m = re.search('Chr(\\d+)\\s+(\\d+)\\s+(\\d+)', line)
                if m:
                    chrom = 'chr{}'.format(m.group(1))
                    start = int(m.group(2)) + 1
                    stop = int(m.group(3))
                    gs_[chrom] = gs_.get(chrom, 0) + stop - start
                    nuc_len += stop - start
                    centromere.append((chrom, start, stop))
            genomesize = gs_
            # for chrom in gs_.keys():
            #     print('{}\t{}\t{}'.format(chrom, gs_[chrom], genomesize[chrom]))
    # filter non-nucleus reads

    nucbams = []
    for fn in set(fns):
        fn_tmp = tempfile.mktemp('.bam')
        cmd = ['samtools', 'view', '-o', fn_tmp, fn]
        for chrm, start, stop in centromere:
            cmd.append('{}:{}-{}'.format(chrm, max(1, start), stop))
        # for ch, ln in genomesize.items():
        #     if re.match('chr\\d$', ch):
        #         cmd.append('{}:1-{}'.format(ch, ln))
        sys.stderr.write(' '.join(cmd) + '\n')
        subprocess.Popen(cmd).wait()
        if os.path.exists(fn_tmp) and os.path.getsize(fn_tmp) > 1000:
            nucbams.append(fn_tmp)
        else:
            sys.stderr.write('failed to write {}\n'.format(fn_tmp))

    cmd = ['macs2', 'callpeak', '-f', 'BAMPE', '--outdir', dstdir, '-g', str(nuc_len), '--shift', '100', '--extsize', '200', '-n', name, '-t'] + nucbams
    sys.stderr.write(' '.join(cmd) + '\n')

    subprocess.Popen(cmd).wait()

    for fn in nucbams:
        if os.path.exists(fn): 
            os.unlink(fn)

        



#    cmd = ['macs2', 'callpeak', 
