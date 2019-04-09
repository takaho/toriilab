import os, sys, re
import pathlib, argparse
import numpy as np
import collections

class MEMEMotif(object):
    def __init__(self, name, length, n_sites=0, score=0, matrix=None):
        self.name = name
        self.length = length
        self.n_sites = n_sites
        self.score = score
        self.set_matrix(matrix)
        self.background = [.25] * 4
        # self.matrix = matrix
    def set_background(self, bg):
        if len(bg) == 4:
            self.background = bg
    def set_matrix(self, mat):
        if mat is None:
            return
        if len(mat) != self.length:
            raise Exception('matrix shape should be [{},4]'.format(self.length))
        for r in mat:
            if len(r) != 4:
                raise Exception('matrix shape should be [{},4]'.format(self.length))
        self.matrix = mat
    def __repr__(self):
        ostr = '{} Score:{:.2e} Sites:{}\n'.format(self.name, self.score, self.n_sites)
        for i in range(4):
            line = ''
            for j in range(self.length):
                if j > 0: line += ' '
                line += '{:.3f}'.format(self.matrix[j][i])
            ostr += line + '\n'
        return ostr
    def output_meme(self, ost):
        a, c, g, t = self.background
        ost.write('MEME version 5\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies (from uniform background):\nA {:.4f} C {:.4f} G {:.4f} T {:.4f}\n\n'.format(a, c, g, t))
        ost.write('MOTIF {}\n\nletter-probability matrix: alength= 4 w= {} nsites= {} E={:.3e}\n\n'.format(self.name, self.length, self.n_sites, self.score))
        for i in range(4):
            line = ''
            for j in range(self.length):
                if j > 0: line += ' '
                line += '{:.3f}'.format(self.matrix[j][i])
            ost.write(line)
            ost.write('\n')
    @classmethod
    def load_dreme(cls, fn):
        if os.path.isdir(fn):
            title = os.path.basename(fn)
            fn = os.path.join(fn, 'dreme.txt')
        elif os.path.isfile(fn):
            title = os.path.basename(os.path.dirname(fn))
        else:
            raise Exception('Skip {} : directory or dreme.txt is acceptable\n'.format(fn))

        if os.path.exists(fn) is False:
            raise Exception('no {}\n'.format(fn))

        with open(fn) as fi:
            contents = fi.readlines()
        num_lines = len(contents)
        i = 0
        motifs = []
        background = [.25, .25, .25, .25]
        while i < num_lines:
            line = contents[i]
            if line.startswith('Background letter frequencies'):
                m = re.match('A\\s+([\\d\\.]+)\\s+C\\s+([\\d\\.]+)\\s+G\\s+([\\d\\.]+)\\s+T\\s+([\\d\\.]+)\\s+', contents[i+1])
                if m:
                    acgt = [float(x_) for x_ in m.groups()]
                    if len(acgt) == 4:
                        background = [float(x_) / sum(acgt) for x_ in acgt]
                i += 1
            elif line.startswith('MOTIF'):
                m = re.match('MOTIF\\s+(\\w+)', line)
                name = m.group(1)
                start = -1
                j = i + 1
                while j < num_lines:
                    if contents[j].startswith('letter-probability'):
                        start = j
                        break
                    j += 1
                if start < 0:
                    i += 1
                    continue
                info = contents[start]
                m = re.search('alength=\\s*(\\d+) w=\\s*(\\d+) nsites=\\s*(\\d+) E=\\s*([\\d\\.e\\-]+)', info) #7.2e-543', info)
                l = int(m.group(1))
                w = int(m.group(2))
                n = int(m.group(3))
                score = float(m.group(4))
                j = start + 1
                matrix = []
                while j < min(start + 1 + w, num_lines):
                    freq = [float(x_) for x_ in re.split('\\s+', contents[j].strip())]
                    matrix.append(freq)
                    j += 1
                motif = MEMEMotif(name, w, n, score, matrix)
                # motif.output_meme(sys.stdout)
                motifs.append(motif)
                i = j + 1
            else:
                i += 1
        return motifs

def check_match(seq1, seq2):
    import Bio.Seq
    f = seq1
    r = Bio.Seq.reverse_complement(seq1)
    if f.find(seq2) >= 0 or r.find(seq2) >= 0 or seq2.find(f) >= 0 or seq2.find(r) >= 0:
        return True
    return False

parser = argparse.ArgumentParser()
parser.add_argument('-i', nargs='+')
parser.add_argument('-o', default='out')

args = parser.parse_args()

# motifs = MEMEMotif.load_dreme('500/GC1.dreme')
# seq = 'TGGGCC'
# for key in motifs:
#     print('{}\t{}\t{}'.format(seq, key.name, check_match(key.name, seq)))
# exit()

    
# TGGGCC
for tf in args.i:
    src_dreme = tf + '.dreme/dreme.txt'
    src_tomtom = tf + '.tomtom/tomtom.tsv'
    src_homer = tf + '.homer/homerMotifs.all.motifs'
    tomtom = collections.OrderedDict()
    homer = collections.OrderedDict()
    dreme = collections.OrderedDict()
    for mot in MEMEMotif.load_dreme(src_dreme):
        dreme[mot.name] = mot

    with open(src_tomtom) as fi:
        for line in fi:
            items = line.strip().split('\t')
            if len(items) < 5: continue
            name = items[0]
            dbname = items[1]
            if name not in tomtom: tomtom[name] = []
            tomtom[name].append(dbname)

    # print(tomtom)
    with open(src_homer) as fi:
        for line in fi:
            if line.startswith('>'):
                items = re.split('\\s+', line[1:-1])
                name = items[0]                
                items = line.strip().split('\t')
                """T:1233.0(15.64%),B:1494.2(6.54%),P:1e-172"""
                prop = {}
                for info in items[5].split(','):
                    m = re.match('(\\w+):([\\d\\.\\-e]+)(\\(([\\d\\.]+)%\\))?', info)
                    if m:
                        field = m.group(1)
                        value = float(m.group(2))
                        if m.group(3) is not None:
                            percentage = float(m.group(4))
                        else:
                            percentage = None
                        prop[field] = value, percentage
                reproduced = '.'
                if 'T' in prop and 'B' in prop:
                    target = prop['T'][1]
                    bg = prop['B'][1]
                    enrichment = np.log2((prop['T'][1] + .1)/(prop['B'][1] + .1))
                else:
                    enrichment = '-'
                pvalue = prop.get('P', [1.0,1.0])[0]
                homer[name] = [name, enrichment, pvalue]
    for name in sorted(homer, key=lambda n:homer[n][1], reverse=True):
        name, enrichment, pvalue = homer[name]
        found = False
        for key, names in tomtom.items():
            if check_match(name, key):
                reproduced = ','.join(names)
                found = True
                break
        if not found:
            # print(key, ','.join(dreme.keys()))
            for key, motif in dreme.items():
                if check_match(name, key):
                    reproduced = 'DREME:{}'.format(mot.name)
                    found = True
                    break
        if found:            
            print('{}\t{}\t{:.2f}\t{:.2e}'.format(name, reproduced, enrichment, pvalue))#str(prop), reproduced))

            



