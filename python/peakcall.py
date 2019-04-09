import os, sys, re, subprocess
import datetime

fn_log = 'log.txt'
touched = {}
data = {}
genomesize = {}

if os.path.exists(fn_log):
    with open(fn_log) as fi:
        for line in fi:
            items = line.split('\t')
            name = items[0]
            data[name] = {}
            for item in items[1:]:
                m = re.match('(.+)?:(\\d+)$', item)
                if m:
                    data[name][m.group(1)] = int(m.group(2))
dstdir = 'autosome'
peakdir = 'macs'
for dn in dstdir, peakdir:
    if os.path.exists(dn) is False:
        os.makedirs(dn)

for fn in sys.argv[1:]:
    if os.path.exists(fn) is False or fn.endswith('.bam') is False: continue
    if os.path.getsize(fn) < 100000000:
        sys.stderr.write('{} has too few reads\n'.format(fn))
        continue
    if re.search('\\.\\d{4}\\.bam$', fn):
        continue
    name = os.path.basename(fn).split('.')[0]
    sys.stderr.write('[{}] processing {}\n'.format(datetime.datetime.now().strftime('%H:%M:%S'), name))
    if len(genomesize) == 0:
        cmd = 'samtools', 'view', '-H', fn
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        for line in proc.stdout:
            m = re.match('@SQ\\s+SN:(\\S+)\\s+LN:(\\d+)', line)
            if m:
                genomesize[m.group(1)] = int(m.group(2))
        #print(genomesize)
        proc.stdout.close()
        proc.wait()
    fn_mac = os.path.join(peakdir, name + '_peaks.xls')
    if os.path.exists(fn_mac) is False:
        # Index
        bai = fn + '.bai'
        if os.path.exists(bai) is False:
            cmd = 'samtools', 'index', fn
            sys.stderr.write('Indexing {}...         \r'.format(name))
            subprocess.Popen(cmd).wait()
            pass
        
        # Generate partial BAM
        nuc = []
        fn_nuc = os.path.join(dstdir, name + '.autosome.bam')
        gs = 0
        if os.path.exists(fn_nuc) is False:
            cmd = ['samtools', 'view', '-o', fn_nuc, fn]
            for chrom, length in genomesize.items():
                if re.match('chr\\d+$', chrom):
                    gs += length
                    cmd.append('{}:1-{}'.format(chrom, length))
            sys.stderr.write(' '.join(cmd) + '\n')
            proc = subprocess.Popen(cmd).wait()

        #peak call
        cmd = 'macs2', 'callpeak', '-f', 'BAMPE', '-g', str(gs), '--outdir', peakdir, '-n', name, '-t', fn_nuc
        sys.stderr.write(' '.join(cmd) + '\n')
        proc = subprocess.Popen(cmd).wait()

        if os.path.exists(fn_nuc):
            os.unlink(fn_nuc)
        
    # exit()

    # cmd = 'samtools', 'view', fn
    # sys.stderr.write(' '.join(cmd) + '\n')
    # proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    # counts = {}
    # for line in proc.stdout:
    #     items = line.decode('ascii').split('\t', 3)
    #     chrom = items[2]
    #     if chrom not in counts:
    #         counts[chrom] = 1
    #     else:
    #         counts[chrom] += 1
    # ostr = name
    # for chrom in sorted(counts.keys()):
    #     ostr += '\t{}:{}'.format(chrom, counts[chrom])
    # with open(fn_log, 'a') as fo:
    #     print(ostr)
    #     fo.write(ostr + '\n')
    # proc.stdout.close()
    # proc.wait()
    
    
"""
            @HD	VN:1.0	SO:coordinate
@SQ	SN:chr1	LN:30427671
@SQ	SN:chr2	LN:19698289
@SQ	SN:chr3	LN:23459830
@SQ	SN:chr4	LN:18585056
@SQ	SN:chr5	LN:26975502
@SQ	SN:chrC	LN:154478
@SQ	SN:chrM	LN:366924
"""
