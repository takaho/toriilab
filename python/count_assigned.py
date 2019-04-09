import os, sys, re



def get_counts(filename):
    """
6:intergenic
0:TSS
1:exon
2:intron
3:5'UTR
4:3'UTR
5:transposon

    """
    counts = [0] * 7
    with open(filename) as fi:
        for line in fi:
            items = line.strip().split('\t')
            t = items[4]
            if t == '.':
                code = 6
            elif t.find(':TSS:') > 0:
                code = 0
            elif t.find(':5\'UTR:') > 0:
                code = 3
            elif t.find(':3\'UTR:') > 0:
                code = 4
            elif t.find(':exon:') > 0:
                code = 1
            elif t.find(':transposon:') >= 0:
                code = 5
            elif t.find(':intron:') >= 0:
                code = 2
            else:
                code = 6
                sys.stderr.write(line)
            counts[code] += 1
    return counts
                
data = {}
for fn in sys.argv[1:]:
    name = os.path.join(fn).split('.')[0]
    #print(name)
    data[name] = get_counts(fn)

labels = {
6:'intergenic',
0:'TSS',
1:'exon',
2:'intron',
3:'5\'UTR',
4:'3\'UTR',
5:'transposon',
}

names = list(data.keys())
print('Group\t' + '\t'.join(names))
for i in range(7):
    ostr = labels[i]
    for name in names:
        ostr += '\t{}'.format(data[name][i])
    print(ostr)

print('Group\t' + '\t'.join(names))
for i in range(7):
    ostr = labels[i]
    for name in names:
        ostr += '\t{:.1f}'.format(data[name][i] * 100.0 / sum(data[name]))
    print(ostr)
#for name in sorted(data.keys()):
    
