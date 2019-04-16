import os, re, sys
import numpy as np
class ATACMotif(object):
    def __init__(self, title, consensus=None):
        self.title = title
        if consensus is None: consensus = title
        self.consensus = consensus
        self.matrix = None
        self.__known = {}
    
    def __give_1(self):
        return 1
    def __give_0(self):
        return 0
    enrichment = property(__give_1)
    pvalue = property(__give_1)
    n_sites = property(__give_0)

    def get_shared_consensus(self, mot):
        shared = []
        for c1 in self.__known.keys():
            if c1 in mot.__known:
                shared.append(c1)
        return shared

    def set_matrix(self, mat):
        if mat is None:
            return
        m = np.array(mat)
        s = m.shape
        if s[0] == 4 and s[1] != 4:
            self.matrix = m.T
        elif s[1] == 4:
            self.matrix = m
        else:
            raise Exception('invalid matrix shape')

    def add_known(self, name, detected_consensus, motif_consensus, score=0):
        self.__known[name] = (detected_consensus, motif_consensus, score)
    
    def get_known(self, namelist=None):
        ostr = ''
        if isinstance(namelist, str):
            namelist = (namelist, )
        for name in sorted(self.__known.keys(), key=lambda k:self.__known[k][1]):
            if namelist is None or name in namelist:
                if ostr != '': ostr += ','
                seq, cns, score = self.__known[name]
                ostr += '{}({})'.format(name, seq)
        if ostr == '': ostr = '.'
        return ostr

    @classmethod
    def __normalize_motif(cls, seq):
        rev = ''
        comp = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'N', '.':'N', 
            'V':'B', 'H':'D', 'D':'H', 'B':'V', 'R':'R', 'W':'W', 'S':'S', 'M':'K', 'K':'M','Y':'R', 'R':'Y', 'U':'A'}
        for i in range(len(seq)):
            n = seq[len(seq) - 1 - i]
            rev += comp.get(n, n)
        return rev if rev < seq else seq
    @classmethod
    def annotate_motifs(cls, motifs, filename, pvalue_threshold=0.05):
        with open(filename) as fi:
            fi.readline()
            for line in fi:
                items = line.strip().split('\t')
                qid = items[0]
                found = False
                for motif in motifs:
                    if motif.title == qid:
                        known = items[1]
                        pvalue = float(items[3])
                        qvalue = float(items[5])
                        evalue = float(items[4])
                        if pvalue < pvalue_threshold:# or 1:
                            motif_consensus = cls.__normalize_motif(items[8])
                            detected = cls.__normalize_motif(items[7])
                            motif.add_known(known, detected, motif_consensus, evalue)
                            found = True
        return motifs
    @classmethod 
    def classify_and_select_best(cls, motifs):
        dbnames = {}
        integrated = {}
        novel_motifs = {}
        sequences = set()

        # aggreagate motifs by binding sequence
        for m in motifs:
            sequences.add(m.consensus)
            novel = True
            for key, val in m.__known.items():
                consensus, motif_consensus, score = val
                if consensus not in integrated:
                    integrated[consensus] = set()
                integrated[consensus].add(key)
                if key not in dbnames:
                    dbnames[key] = []
                dbnames[key].append(m)
                novel = False
            if novel:
                dbnames[m.title] = [m, ]

        most_hits = {}
        touched = set()
        for name, mots in dbnames.items():
            joined = None
            for mot in mots:
                if mot.consensus in touched: continue
                for con, names in integrated.items():
                    # if con in touched: continue
                    if name in names:
                        joined = '//'.join(sorted(names))
                        touched.add(mot.consensus)
                        break
            if joined is not None:
                most_hits[joined] = list(sorted(dbnames[name], key=lambda m:m.n_sites, reverse=True))[0]
            # else:
            #     most_hits[name] = list(sorted(mots, key=lambda m:m.n_sites, reverse=True))[0]
                # print('new mot', name, most_hits[name].consensus)

        return most_hits

class MEMEMotif(ATACMotif):
    def __init__(self, name, consensus, length, n_sites=0, score=0, matrix=None):
        ATACMotif.__init__(self, name, consensus)
        self.name = name
        self.length = length
        self.__hits = n_sites
        self.score = score
        self.set_matrix(matrix)
        self.background = [.25] * 4
        # self.matrix = matrix
    def set_background(self, bg):
        if len(bg) == 4:
            self.background = bg
    def __get_hits(self):
        return [self.n_sites, 1.0]
    n_sites = property(lambda s:s.__hits)
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
    def load_motifs(cls, fn):
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
                m = re.match('MOTIF\\s+(\\w+)\\s+(\\S+)', line)
                name = m.group(2)
                consensus = m.group(1)
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
                # print(name, consensus)
                motif = MEMEMotif(consensus, consensus, w, n, score, matrix)
                motifs.append(motif)
                i = j + 1
            else:
                i += 1
        return motifs

class HOMERMotif(ATACMotif):
    def __init__(self, motif, name, matrix, info):
        ATACMotif.__init__(self, motif, motif)
        self.motif = motif
        self.name = name
        self.dbname = None
        # self.matrix = matrix
        self.info = info
        self.set_matrix(matrix)
    def set_knwon_name(self, name):
        if self.dbname is None:
            self.dbname = []
        self.dbname.append(name)

    def __get_enrichment(self):
        if 'T' in self.info and 'B' in self.info:
            t = self.info['T'][1]
            b = self.info['B'][1]
            if t > 0 and b > 0:
                return t / b
        return 1.0

    def __get_hits(self):
        if 'T' in self.info:
            return self.info['T'][0]
        else:
            return 0

    def __get_background(self):
        if 'B' in self.info:
            return self.info['B'][0]
        else:
            return 0

    def __get_known_name(self):
        priority_names = 'SEP1', 'SOC1', 'PI', 'AGL15', 'ATHB-5', 'SEP3', 'TBP3', 'SEP4', 'SEP5',  
        if self.dbname is None:
            return '.'
        for n in priority_names:
            if n in self.dbname:
                return n
        return self.dbname[0]

    enrichment = property(__get_enrichment)
    hits = property(__get_hits)
    background = property(__get_background)
    known = property(__get_known_name)
    n_sites = property(__get_hits)


    def __repr__(self):
        ostr = '{} {} {} E:{},T:{},B:{}\n'.format(self.motif, self.name, self.known, self.enrichment, self.hits, self.background)
        for i in range(4):
            line = ''
            for j in range(self.length):
                if j > 0: line += ' '
                line += '{:.3f}'.format(self.matrix[j][i])
            ostr += line + '\n'
        return ostr
    @classmethod
    def convert_to_meme(cls, motifs, filename_template=None):
        if filename_template is None:
            filename_template = os.path.join(os.path.dirname(__file__), 'meme_template.txt')
        with open(filename_template) as fi:
            contents = fi.read()
        ostr = contents
        for i, motif in enumerate(motifs):
            ostr += 'MOTIF {} HOMER-{}\n\n'.format(motif.motif, i + 1)#motif.name)
            ostr += 'letter-probability matrix: alength= 4 w= {} nsites= {} E= {:.4e}\n'.format(motif.length, motif.info.get('T', [0,0])[0], motif.info.get('P', 0))
            for row in motif.matrix:
                ostr += '{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\n'.format(row[0], row[1], row[2], row[3])
            ostr += '\n\n'
        return ostr

    length = property(lambda s:len(s.matrix))
    @classmethod
    def load_motifs(cls, srcdir):
        motifs = {}
        for fn in os.listdir(srcdir):
            m = re.match('motif(\\d+)\\.motif$', fn)
            if m is None:
                continue
            number = int(m.group(1))
            with open(os.path.join(srcdir, fn)) as fi:
                # items = re.split('\\s+', fi.readline()[1:-1])
                items = fi.readline()[1:-1].split('\t')
                # print(items)
                motif = items[0]
                name = items[1]
                evalue = float(items[2])
                score = float(items[3])
                dig = float(items[4])
                info = {}
                for elem in items[5].split(','):
                    key, val = elem.split(':', 1)
                    if key == 'T' or key == 'B':
                        m = re.match('([\\d\\.]+)\\(([\\d\\.]+)%\\)', val)
                        if m:
                            count, percentage = float(m.group(1)), float(m.group(2))
                            info[key] = count, percentage
                    elif key == 'P':
                        pvalue = float(val)
                info['E'] = evalue
                matrix = []
                for line in fi:
                    items = re.split('\\s+', line.strip())
                    matrix.append([float(x) for x in items])
                motifs[number] = HOMERMotif(motif, name, matrix, info)
        return [motifs[n] for n in sorted(motifs.keys())]

def load_motifs(dn):
    if os.path.isdir(dn) and os.path.exists(os.path.join(dn, 'motif1.motif')):
        motifs = HOMERMotif.load_motifs(dn)
        tomtom = os.path.join(os.path.dirname(os.path.dirname(dn)), 'tomtom', 'tomtom.txt')
    elif os.path.isdir(dn) and os.path.exists(os.path.join(dn, 'dreme.txt')):
        motifs = MEMEMotif.load_motifs(os.path.join(dn, 'dreme.txt'))
        tomtom = os.path.join(dn, 'tomtom', 'tomtom.txt')
    elif os.path.basename(dn) == 'dreme.txt':
        motifs = MEMEMotif.load_motifs(dn)   
        tomtom = os.path.join(os.path.dirname(dn), 'tomtom', 'tomtom.txt')
    else:
        sys.stderr.write('cannot process {}'.format(dn))
        raise Exception("")
    if os.path.exists(tomtom):
        ATACMotif.annotate_motifs(motifs, tomtom)
    else:
        sys.stderr.write('no tomtom {}\n'.format(tomtom))    
    return motifs

if __name__ == '__main__':
    #import argparse
    #parser = argparse.ArgumentParser()
    dn = sys.argv[1]
    if os.path.isdir(dn) and os.path.exists(os.path.join(dn, 'motif1.motif')):
        motifs = HOMERMotif.load_motifs(dn)
        tomtom = os.path.join(os.path.dirname(os.path.dirname(dn)), 'tomtom', 'tomtom.txt')
    elif os.path.isdir(dn) and os.path.exists(os.path.join(dn, 'dreme.txt')):
        motifs = MEMEMotif.load_motifs(os.path.join(dn, 'dreme.txt'))
        tomtom = os.path.join(dn, 'tomtom', 'tomtom.txt')
    elif os.path.basename(dn) == 'dreme.txt':
        motifs = MEMEMotif.load_motifs(dn)   
        tomtom = os.path.join(os.path.dirname(dn), 'tomtom', 'tomtom.txt')
    else:
        sys.stderr.write('cannot process {}'.format(dn))
        raise Exception("")
    if os.path.exists(tomtom):
        ATACMotif.annotate_motifs(motifs, tomtom)
    else:
        sys.stderr.write('no tomtom {}\n'.format(tomtom))    

    arranged = ATACMotif.classify_and_select_best(motifs)
    for name, m in arranged.items():
        print('{}\t{}\t{:.2f}\t{:.0f}'.format(name, m.consensus, m.enrichment, m.n_sites))

def load_consensus(filename):
    def convert_freq_to_base(freq):
        index = np.argmax(freq)
        if freq[index] < 0.30:
            return 'N'
        elif freq[index] > .9:
            return 'ACGT'[index]
        #mindex = np.argmin(freq)
        # vals = list(sorted(freq))
        flag = 0
        for i, v in enumerate(freq):
            if v > .3:
                flag |= 1 << i
        return 'NACMGRSTWYHKDBN'[flag]
        
    motifs = {}
    with open(filename) as fi:
        name = matrix = None
        while 1:
            line = fi.readline()
            if line == '': break
            if line.startswith('MOTIF'):
                name = re.split('\\s+', line.strip())[-1]
            elif name is not None and line.startswith('letter-probability matrix:'):
                info = line.split(':', 1)[1]
                prop = {}
                m = re.search('w=\\s*(\\d+)', info)
                if m:
                    length = int(m.group(1))
                else:
                    length = 0
                seq = ''
                for i in range(length):
                    items = re.split('\\s+', fi.readline().strip())
                    if len(items) < 4:
                        continue
                    seq += convert_freq_to_base([float(x_) for x_ in items[0:4]])
                if len(seq) > 3:
                    motifs[name] = seq
                name = None
    return motifs
