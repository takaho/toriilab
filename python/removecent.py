import sys
import re

regions = {}
with  open('TAIR10_nocentromere.bed') as fi:
    for line in fi:
        m = re.search('Chr(\\d+)\\s+(\\d+)\\s+(\\d+)', line)
        if m:
            chrom = 'chr' + m.group(1)
            if chrom not in regions: regions[chrom] = []
            regions[chrom].append((int(m.group(2)), int(m.group(3))))

span_mask = 0
genome = {}
with open('/mnt/smb/tae/tair10/tair10.fa') as fi:
    chrom = seq = ''
    for line in fi:
        if line.startswith('>'):
            if seq != '': genome[chrom] = seq
            chrom = line[1:-1].strip()
            seq = ''
        else:
            seq += line.strip()
    if seq != '': genome[chrom] = seq
for chrom in regions.keys():
    r0, r1 = regions[chrom]
    seq = genome[chrom]
    span_mask += (r1[0] - r0[1])
    sys.stderr.write('{} {} {}\n'.format(r0, r1, span_mask))
    masked = seq[r0[0]:r0[1]] + 'N' * (r1[0] - r0[1]) + seq[r1[0]:r1[1]]
    print('>{}\n{}'.format(chrom, masked))
    sys.stderr.write('{}\t{}\n'.format(chrom, span_mask))
        

            
            
