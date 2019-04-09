import os, sys, re
# import reportlab.pdfgen.canvas
import pathlib

class ATACMotif(object):
    def __init__(self, name, sequence, pwm, pvalue=1.0, score=0, enrichment=1.0):
        self.name = name
        self.sequence = sequence
        self.pvalue = pvalue
        self.enrichment = enrichment
        self.score = score
        self.pwm = pwm
    def __repr__(self):
        return '{}\t{}\tn:{}\tP:{:.2e}\tE:{:.2f}'.format(self.sequence, self.name, len(self.pwm), self.pvalue, self.enrichment)
    @classmethod
    def load_dreme(cls, filename):
        motif = ''
        pwm = None
        score = 0
        pval = 1.0
        motifs = []
        with open(filename) as fi:
            for line in fi:
                # line = line.decode('ascii')
                m = re.match('MOTIF\s+(\\w+)\\sDREME', line)
                if m:
                    if pwm is not None and len(pwm) > 0:
                        motifs.append(ATACMotif(motif, motif, pwm, pval, score))#(motif, score, pval, pwm))
                    name = motif = m.group(1)
                    score = 0
                    pwm = None
                    continue
                if line.startswith('# BEST'):
                    items = re.split('\\s+', line.strip())
                    pval = float(items[-2])
                    score = float(items[-1])
                if line.startswith('letter-probability matrix:'):
                    # m = re.search('E=\\s+(.*)')
                    # if m:
                    #     score = float(m.group(1))
                    pwm = []
                elif pwm is not None:
                    items = re.split('\\s', line.strip())
                    if len(items) == 4:
                        pwm.append([float(x) for x in items])
        if pwm is not None and len(pwm) > 0:
            motifs.append(ATACMotif(name, motif, pwm, pval, score))#(motif, score, pval, pwm))
        return motifs
    @classmethod
    def load_homer(cls, dirname):
        motifs = []
        score = 0
        pval = 1.0
        for fn in os.listdir(dirname):
            m = re.match('motif(\\d+)\\.motif$', fn)
            if m is None: continue
            num = int(m.group(1))
            with open(os.path.join(dirname, fn)) as fi:
                header = next(fi)
                items = re.split('\\s+', header[1:-1])
                name = re.sub('/Homer\\(.*', '', items[1].split(',')[1].split(':')[1])
                #.T:72.0(15.00%),B:267.1(0.58%),P:1e-75
                motif = items[0]
                for elem in items[-1].split(','):
                    mv = re.match('(\\w+):([e\\-\\d\\.]+)(\\(([\\d\\.]+)%\\))?', elem)
                    if mv:
                        # print(mv.groups())
                        key = mv.group(1)
                        value = float(mv.group(2))
                        if key == 'P':
                            pval = value
                        elif key == 'T':
                            target = float(mv.group(4))
                        elif key == 'B':
                            background = float(mv.group(4))
                enrichment = target / background if background > 0 else target
                pwm = []
                for line in fi:
                    pwm.append([float(x_) for x_ in line.strip().split('\t')])
                motifs.append(ATACMotif(name, motif, pwm, pval, enrichment, enrichment))

        return motifs
    @classmethod
    def load_motifs(cls, srcdir):
        path = pathlib.Path(srcdir)
        motifs = None
        if path.name == 'dreme.txt':
            motifs = cls.load_dreme(fitem.as_posix())
        if path.name == 'homerResults':
            motifs = cls.load_homer(fitem.as_posix())
        else:
            for fitem in path.iterdir():
                if fitem.name == 'dreme.txt':
                    motifs = cls.load_dreme(fitem.as_posix())
                if fitem.name == 'homerResults':
                    motifs = cls.load_homer(fitem.as_posix())
        if motifs is None:
            raise Exception('no motifs included')
        return sorted(motifs, key=lambda m:m.pvalue)

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', default='peakcall')
parser.add_argument('-n', nargs='+', default="AtML1,SPCH,MUTE,FAMA,GC1")
parser.add_argument('--graph', default='motif.pdf')
parser.add_argument('-o', default=None)

args = parser.parse_args()


names = {}
results = []
experiments = []
name2seq = {}
basedir = args.i
import motif
tfs = args.n.split(',')
#tfs = 'AtML1','SPCH','MUTE','FAMA','GC1'
for tf in tfs:
    homer = ATACMotif.load_motifs(os.path.join(basedir, 'homer_{}'.format(tf)))#peakcall/homer_{}'.format(tf))#fn)
    meme = ATACMotif.load_motifs(os.path.join(basedir, 'meme_{}'.format(tf)))#peakcall/dreme_{}'.format(tf))#fn)
    mememotifs = []
    homermotifs = []
    for m in meme:
        p = motif.PWM(m.name, m.sequence, 'Arabidopsis')
        mat = [[], [], [], []]
        for row in m.pwm:
            for i in range(4):
                mat[i].append(row[i])
        p.set_matrix(mat)
        mememotifs.append(p)
    experiments.append(tf)
    detected = {}
    for hm in homer:
        accepted = False
        if hm.enrichment < 2: continue
        for m in mememotifs:
            score = m.evaluate_weight_sequence(hm.pwm)
            if score > 0.6:
                accepted = True
                print('Accepted by {} : {} {}'.format(m.name, score, str(hm)))
                homermotifs.append(hm)
                name2seq[hm.name] = hm.sequence
                detected[hm.name] = hm
                names[hm.name] = names.get(hm.name, 0) + 1
                break
        if not accepted:
            print('{} was rejected'.format(hm.name))
    results.append(detected)
#     break

# exit()
    
    
# for fn in sys.argv[1:]:
#     motifs = ATACMotif.load_motifs(fn)
#     print(fn)
#     experiments.append(fn)
#     # results.append(motifs)
#     detected = {}
#     for i, m in enumerate(motifs):
#         if m.enrichment < 2: continue
#         name2seq[m.name] = m.sequence
#         detected[m.name] = m
#         print('{}:{}'.format(i + 1, str(m)))
#         names[m.name] = names.get(m.name, 0) + 1
#     results.append(detected)
emat = {}
accum_pwm = {}
for name in sorted(names.keys(), key=lambda n:names[n], reverse=True):
#    ostr = name
    ostr = name2seq[name]
    row = []
    for result in results:
        if name in result:
            e = result[name].enrichment
            ostr += '\t{:.2f}'.format(result[name].enrichment)
            if name not in accum_pwm:
                accum_pwm[name] = []
            accum_pwm[name].append(result[name].pwm)
        else:
            e = 1.0
            ostr += '\t.'
        row.append(e)
    emat[name] = row, ostr
#    print(ostr)
def get_order(enrichment):
    score = 0
    for i, e in enumerate(enrichment):
        score *= 10
        if e > 1:
            score += min(4, e)
    return score

import numpy as np
import reportlab.pdfgen.canvas

cnv = reportlab.pdfgen.canvas.Canvas(args.graph)

left, top = 100, 700
motif_height = 12
motif_index = 0
balloon_left = left + 150
balloon_width = 20
cnv.setFont('Helvetica', 8)
for i, label in enumerate(tfs):
    cnv.drawString(balloon_left + i * balloon_width - 10, top + 20, label)
    pass
#
print('Sequence\t' + '\t'.join(tfs) + '\tName')


def get_enrichment_radius(enrichment):
    if enrichment > 1:
        return np.log2(enrichment) * 2
    return 1.0

for name in sorted(emat.keys(), key=lambda n:get_order(emat[n][0]), reverse=True):
    print('{}\t{}'.format(emat[name][1], name))
    pwms = accum_pwm[name]
    if len(pwms) == 1 or 1:
        matrix = np.array(pwms[0]).T
    # forward or complementary?
    fm = rm = 0
    for i, base in enumerate(name2seq[name]):
        if base == 'A':
            index = 0
        elif base == 'C':
            index = 1
        elif base == 'G':
            index = 2
        elif base == 'T':
            index = 3
        else:
            continue
        rev = 3 - index
        fm += matrix[index][i]
        rm += matrix[rev][matrix.shape[1] - 1 - i]
    if fm < rm:
        matrix = np.flip(np.flip(matrix, 1), 0)
    y = top - (motif_height + 5) * motif_index
    freq = matrix.T
    # print(freq)
    # print(matrix.T)
    motif.SequenceLogo.draw_on_pdf(canvas=cnv, frequencies=matrix.T, position=(left, y), size=(matrix.shape[1] * 10, motif_height))
    items = emat[name][1].split('\t')
    seq = items[0]
    for i in range(5):
        if items[i + 1] == '.':
            enr = 0
        else:
            enr = float(items[i + 1])
        cnv.saveState()
        cnv.translate(balloon_left + i * balloon_width, y)
        if enr > 1:
            radius = get_enrichment_radius(enr)#min(enr, balloon_width / 2)
            cnv.scale(0.666, 1)#, y, radius, 1, 1)
            cnv.setFillColor('cadetblue')
            cnv.circle(0, 0, radius, 1, 1)
        else:
            cnv.setFillColor('gray')
            cnv.circle(0, 0, 1, 0, 1)
        #cnv.circle(balloon_left + i * balloon_width, y, radius, 1, 1)
        cnv.restoreState()
    cnv.setFillColor('black')
    cnv.drawString(balloon_left + 6 * balloon_width, y, name)
    # print(np.array(matrix * 100, dtype=np.int))
    motif_index += 1

enr = 2
x = balloon_left
while 1:
    radius = get_enrichment_radius(enr)
    if radius > balloon_width / 2:
        break
    cnv.saveState()
    cnv.scale(0.666, 1)#, y, radius, 1, 1)
    cnv.setFillColor('cadetblue')
    cnv.circle(x, 100, radius, 1, 1)
    cnv.restoreState()
    x += balloon_width

    enr *= 2


cnv.save()