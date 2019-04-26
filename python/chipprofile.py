#!/usr/bin/python
#-*-coding:utf-8-*-

"""ChIP-chipのTSS周りの状況を書くためのスクリプトを何度も作成するのが嫌なのでまとめる
"""

import math, sys, types, os, tempfile
#import mousegene
import reportlab.pdfgen.canvas

class Profile(object):
    """Profiles of ChIP-chip"""
    MARGIN_5 = 6000
    MARGIN_3 = 6000
    MAXIMUM_VALUE  = 2 ** 5
    def __init__(self, identifier, symbol=None, alignment=None, status=None ):
        gene = None
        if type(identifier) is types.IntType:
            geneid = identifier
        else:
            gene = mousegene.get_gene_by_name(identifier)
            geneid = gene.geneid
            pass
        self.__geneid = geneid
        if symbol is None:
            if gene is None: gene = mousegene.get_gene_by_geneid(geneid)
            symbol = gene.symbol
            pass
        self.__symbol = symbol
        self.__alignment = alignment
        self.__refseq = None
        self.__status = status
#         if alignment is None:
#             if gene is None: gene = mousegene.get_gene_by_geneid(geneid)
#             self.__alignment = gene
        self.__data = {}
        pass

    def get_maximum(self, left_margin=-4000, right_margin=4000, require_support=False):
        max_pos = None
        max_value = None
        for pos, value in self.__data.items():
            if pos < left_margin or pos > right_margin: continue
            if value is None or value > max_value:
                max_pos = pos
                max_value = value
                pass
            pass
        return max_pos, max_value

    def get_geometric_mean(self, left_margin=-4000, right_margin=4000):
        n = 0
        gm = 0.0
        for pos, value in self.__data.items():
            if pos < left_margin or pos > right_margin: continue
            gm += value
            n += 1
            pass
        #if n <= 1: return None
        return n, gm / n

    def set_data(self, position, logratio, error=None, base=None):
        """Position and log ratio of IP and input.
        This method dicards probes out of MARGINs.
        """
        if position + self.MARGIN_5 < 0 or position > self.MARGIN_3: return
        if error is not None:
            if base is not None:
                logratio = math.log10(math.pow(base, logratio))
                pass
            self.__data[position] = logratio, error
        else:
            if base is not None:
                logratio = math.log10(math.pow(base, logratio))
                pass
            self.__data[position] = float(logratio)
            pass
        pass

    def __get_alignment(self):
        if self.__alignment is None:
            try:
                self.__alignment = mousegene.get_gene_by_geneid(self.__geneid).get_alignment()
            except:
                return "not_aligned"
            pass
        if self.__alignment is not None:
            chrm, ori, start, end = self.__alignment
        else:
            chrm, ori, start, end = '', '', 0, 0
        return "chr%s:%s:%d:%d" % (chrm, ori, start, end)

    def __get_data(self):
        for pos in sorted(self.__data.keys()):
            v = self.__data[pos]
            if type(v) is types.FloatType:
                yield pos, v
            else:
                yield pos, v[0]
            pass
        pass

    def __getitem__(self, pos):
        #if pos in self.__data:
        return self.__data[pos]

    def __contains__(self, pos):
        return pos in self.__data

    def __set_refseq(self, value):
        self.__refseq = value
        pass
    alignment = property(__get_alignment)
    geneid = property(lambda s: s.__geneid)
    symbol = property(lambda s: s.__symbol)
    data = property(__get_data, doc="enumerate position and value for each point")
    refseq = property(lambda s: s.__refseq, __set_refseq)
    status = property(lambda s: s.__status if s.__status is not None else 'protein-coding')
    @classmethod
    def set_range(cls, left, right):
        """Set acceptable range"""
        cls.MARGIN_5 = abs(left)
        cls.MARGIN_3 = abs(right)
        pass

    @staticmethod
    def parzen(dist, size):
        #print(dist, size)
        x = abs(float(dist)) * 2 / size
        if x > 2.0: return 0.0
        if x > 1.0:
            return (2.0 - x) ** 3 * 0.36363867770067627703
        else:
            return (1.0 - 1.5 * x ** 2 + 0.75 * x ** 3) * 1.45455471080270510812;
        pass
    def draw_line(self, canvas, rectangle, valuerange=None, font='Courier', fontsize=9, color='crimson', fill=False, smoothlevel=0, frame=True):
        canvas.saveState()
        left, top, right, bottom = rectangle
        try:
            if valuerange is None:
                valuerange = [[-6000,4000], [0, int(max(self.__data.values()) + 1)]]
            #
            xmin = min(valuerange[0])
            xmax = max(valuerange[0])
            ymin = min(valuerange[1])
            ymax = max(valuerange[1])
            canvas.setFillColor(color)
            xcnv = lambda x: float(x - xmin) / (xmax - xmin) * (right - left) + left
            ycnv = lambda y: float(y - ymin) / (ymax - ymin) * (top - bottom) + bottom
            path = None
            if smoothlevel == 0:
                for p in sorted(self.__data.keys()):
                    if p < xmin: continue
                    if p > xmax: break
                    vx = xcnv(p)
                    if vx > right: break
                    vy = ycnv(self.__data[p])
                    if path is None:
                        path = canvas.beginPath()
                        path.moveTo(vx, vy)#bottom)
                    else:
                        path.lineTo(vx, vy)
            #
            else:
                weights = []#0] * smoothlevel * 2
                for i in range(smoothlevel * 2):
                    weights.append(parzen(abs(i - smoothlevel), smoothlevel))
                #
                dmin = 100000
                prev = None
                positions = sorted(self.__data.keys())
                for p in positions:
                    if prev is not None:
                        dmin = min(dmin, p - prev)
                    prev = p
                #print(dmin)
                #
                p = positions[0]
                while p <= positions[-1]:
                    if p > xmax: break
                    if p >= xmin:
                        w = v = 0.0
                        for i in range(smoothlevel * 2):
                            p_ = p + (i - smoothlevel) * dmin
                            w += weights[i]
                            v += weights[i] * self.__data.get(p_, 0)
                            #print(p_, weights[i], self.__data.get(p_, 0), p_ in self.__data)
                        #
                        vx = xcnv(p)# + dmin)# * 0.5)
                        if vx > right: break
                        v = (v / w) if w > 0 else 0
                        vy = ycnv(min(v, ymax))
                        #print(p, v, w)
                        if path is None:
                            path = canvas.beginPath()
                            path.moveTo(vx, ycnv(ymin))
                        else:
                            path.lineTo(vx, vy)
                    #
                    p += dmin
            #
            if path:
                canvas.setStrokeColor(color)
                canvas.setFillColor(color)
                canvas.drawPath(path, 1, fill)
            #
            if frame:
                canvas.setStrokeColor('black')
                ticks = canvas.beginPath()
                x0 = xcnv(valuerange[0][0])
                y0 = ycnv(valuerange[1][0])
                x1 = xcnv(valuerange[0][1])
                y1 = ycnv(valuerange[1][1])
                ticks.moveTo(x0, bottom)
                x = min(valuerange[0])
                while x <= valuerange[0][1] + 1e-3:
                    ticks.moveTo(xcnv(x), bottom)
                    ticks.lineTo(xcnv(x), bottom - 2)
                    x += 1000
                #
                y = min(valuerange[1])
                if ymax - ymin <= 1: ydiv = 0.25
                elif ymin - ymax > 100: ydiv = 100
                elif ymax - ymin > 10: ydiv = 5
                else: ydiv = 1
                while y <= valuerange[1][1] + 0.0001:
                    ticks.moveTo(left, ycnv(y))
                    ticks.lineTo(left - 2, ycnv(y))
                    y += ydiv
                #
                canvas.drawPath(ticks)
                canvas.rect(left, bottom, right - left, top - bottom)
                canvas.setFillColor('black')
                canvas.setFont(font, fontsize)
                canvas.drawString(left - canvas.stringWidth(repr(ymax)) - 4, top, repr(ymax))
                canvas.drawString(left - canvas.stringWidth(repr(ymin)) - 4, bottom, repr(ymin))
                canvas.setFillColorRGB(.01, .005, .012)
                canvas.drawString(left - canvas.stringWidth(repr(xmin)) * .5, bottom - 4 - fontsize, repr(xmin))
                canvas.drawString(right - canvas.stringWidth(repr(xmax)) * .5 , bottom - 4 - fontsize, repr(xmax))
                canvas.drawString(xcnv(0) - .5 * canvas.stringWidth('0') , bottom - 4 - fontsize, '0')

        except:
            #
            raise
        finally:
            canvas.restoreState()
    def draw(self, canvas, rectangle, font='Courier', fontsize=9, bar_width=4, bar_color='orange', error_color='gray', window_size=0, cpginfo=None):
        """Drag profile on a PDF canvas provided by reportlab library.
        cpginfo : list of (start,end)
        """
        canvas.saveState()
        canvas.setFillColor('black')
        if fontsize > 0:
            canvas.setFont(font, fontsize)
            pass
        left, top, right, bottom = rectangle
        if cpginfo is None:
            info = self.__symbol
            pass
        else:
            info = "GeneID:%d %s %s" % (self.__geneid, self.__symbol, self.alignment)
            pass
        #info = "%s" % (self.__geneid, self.__symbol, self.alignment)
        if fontsize > 0:
            canvas.drawString( left, top, info )
            #canvas.drawString( left, top, self.__symbol)#info )
            #canvas.drawString( left, top -10, self.alignment)#self.__symbol)#info )
            pass
        celing = top - 10
        baseline = bottom + 15
        div = (celing - baseline) / math.log10(self.MAXIMUM_VALUE)
        xconv = lambda x_: (x_ + self.MARGIN_5) / float(self.MARGIN_3 + self.MARGIN_5) * (right - left) + left
        #yconv = lambda y_: baseline if y_ <= 0.0 else baseline + y_ * div
        yconv = lambda y_: baseline + y_ * div
        canvas.setFillColor(bar_color)
        canvas.setStrokeColor(error_color)
        canvas.setLineWidth(0.5)
        #bar_width = 4
        error_width = 1.5
        if window_size <= 0:
            for pos, datum in self.__data.items():
                if type(datum) is types.FloatType:
                    value = datum
                    lower = upper = None
                else:
                    value = datum[0]
                    sd = datum[1]
                    lower = value - sd
                    upper = value + sd
                    pass
                x = xconv(pos)
                y0 = yconv(value)
                if value > 0.0:
                    canvas.setFillColor(bar_color)
                    canvas.rect(x - bar_width / 2, baseline, bar_width, y0 - baseline, 0, 1)
                else:
                    canvas.setFillColor('cyan')
                    #print(x, y0, self.__symbol)
                    canvas.rect(x - bar_width / 2, y0, bar_width, (baseline - y0) , 0, 1)
                    #canvas.rect(x - bar_width / 2, y0, bar_width, (baseline - y0) , 0, 1)
                    pass
                if lower is not None and upper is not None:
                    y1 = yconv(lower)
                    y2 = yconv(upper)
                    canvas.line(x - error_width, y1, x + error_width, y1)
                    canvas.line(x - error_width, y2, x + error_width, y2)
                    canvas.line(x, y1, x, y2)
                    pass
                pass
#                pass
            pass
        else:
            positions = sorted(self.__data.keys())
            path = None
            def window_function(x):
                x = abs(x)
                if x > 1.0: return 0.0
                x *= 2
                if x > 1.0:
                    return math.pow( 2.0 - x, 3 ) * 0.25
                else:
                    return 1.0 - 1.5 * x * x + 0.75 * x * x * x
                pass

            x = -Profile.MARGIN_5
            while x < Profile.MARGIN_3:
                sum_v = sum_w = 0
                for p in positions:
                    dist = abs(x - p)
                    if dist + window_size >= 0 and dist <= window_size:
                        w = window_function(float(dist) / window_size)
                        datum = self.__data[p]
                        value = w * datum if type(datum) is types.FloatType else datum[0]
                        sum_w += w
                        sum_v += value * w
                        pass
                    pass
                x1 = xconv(x)
                #if sum_w == 0:
                y1 = yconv(sum_v/sum_w) if sum_w > 0 and sum_v > 0 else yconv(0.0)
                if path is None:
                    path = canvas.beginPath()
                    path.moveTo(x1,y1)
                else:
                    path.lineTo(x1, y1)
                    pass
                pass
                x += window_size * 0.25
                pass
            if path is not None:
                canvas.setStrokeColor(bar_color)
                canvas.drawPath(path)
                pass
            pass
        y = 1
        canvas.setLineWidth(1.0)
        canvas.setStrokeColor('black')
        canvas.line(left, baseline, left, yconv(math.log10(self.MAXIMUM_VALUE)))
        if self.MAXIMUM_VALUE < 63356:
            while y <= self.MAXIMUM_VALUE:
                y0 = yconv(math.log10(y))
                canvas.line(left, y0, left - 2, y0)
                y *= 2
                pass
            pass
        if cpginfo is None:
            canvas.line(left, baseline, right, baseline)
        else:
            canvas.saveState()
            canvas.setFillColorRGB(1.0, 0.8, 0.8)
            canvas.rect(left, baseline, right - left, 2, 0, 1)#baseline)
            canvas.setFillColor('red')
            for cpgstart, cpgend in cpginfo:
                xc0 = xconv(cpgstart)
                xc1 = xconv(cpgend)
                canvas.rect(xc0, baseline, xc1 - xc0, 2, 0, 1)
                pass
            canvas.restoreState()
            pass
        x = - self.MARGIN_5
        if self.MARGIN_3 - self.MARGIN_5 > 500000:
            tstep = 100000
        elif self.MARGIN_3 - self.MARGIN_5 > 50000:
            tstep = 10000
        elif self.MARGIN_3 - self.MARGIN_5 > 5000:
            tstep = 10000
        else:
            tstep = 1000
            pass
        while x <= self.MARGIN_3:
            x0 = xconv(x)
            canvas.line(x0, baseline, x0, baseline - 2)
            #if x == 0:
            x += tstep
            pass

        x0 = xconv(0)
        canvas.setStrokeColor('gray')
        path = canvas.beginPath()
        path.moveTo(x0, baseline)
        path.lineTo(x0, baseline - 5)
        path.lineTo(x0 + 20, baseline - 5)
        canvas.drawPath(path)
        path = canvas.beginPath()
        path.moveTo(x0 + 20, baseline - 5)
        path.lineTo(x0 + 18, baseline - 3)
        path.lineTo(x0 + 18, baseline - 7)
        path.lineTo(x0 + 20, baseline - 5)
        canvas.setFillColor('gray')
        canvas.drawPath(path, 0, 1)
        #canvas.line(x0, baseline, x0, celing)

        canvas.restoreState()
        pass

    @classmethod
    def load(cls, filename, accepted=None, linear=False):
        """Load ChIP-chip data. If you do not need all genes, designate acceptable genes (in GeneID integer) as the parameter accepted. """
        profiles = {}
        profile = None
        #print(accepted)
        with open(filename) as fi:
            for line in fi:
                if line.startswith('>'):
                    items = line[1:-1].split('\t')
                    if items[0].isdigit():
                        geneid = int(items[0])
                        alignemnt = items[1]
                        symbol = items[2]
                        refseq = '.'
                        status = '.'
                    else:
                        geneid = int(items[0][7:])
                        refseq = items[1]
                        symbol = items[2]
                        alignment = items[3]
                        status = 'protein-coding' if len(items) > 4 else items[4]
                        pass

    #                print(geneid, geneid in accepted)
    #                geneid = int(items[0][7:]) if items[0].isdigit() is False else int(items[0])
                    if (accepted is None or geneid in accepted) and geneid not in profiles:
                        #print(geneid, symbol)
    #                    refseq = items[1]
    #                    symbol = items[2]
    #                    alignment = items[3]
                        try:
                            chrm, ori, start, end = alignment.split(':')
                        except:
                            chrm, ori, start, end = '?', '?', -1,-1
                            pass
                        #chrm, ori, start, end = alignment.split(':')
                        profile = Profile(geneid, symbol, (chrm, ori, int(start), int(end)), status)
                        profile.refseq = refseq
                        profiles[geneid] = profile
                    else:
                        profile = None
                    pass
                elif line.startswith('/'):
                    profile = None
                    pass
                elif profile is not None:
                    items = line.strip().split('\t')
                    if items[0].startswith('A'): items = items[1:]
                    position = int(items[0])
                    value = float(items[2])
                    error = None if (len(items) < 4 or len(items[3]) < 2) else float(items[3])
                    if linear:
                        if value > 0:
                            profile.set_data(position, math.log(value) / math.log(2.0), error)
                            pass
                        pass
                    else:
                        profile.set_data(position, value, error)
                        pass
                    pass
                pass
            pass
        return profiles
    @classmethod
    def read_profile(cls, filename, accepted=None):
        """Load ChIP-chip data. If you do not need all genes, designate acceptable genes (in GeneID integer) as the parameter accepted.
        This method yields each profile on begin loaded.
        """
        #profiles = {}
        profile = None
        touched = set([])
        for line in open(filename):
            if line.startswith('>'):
                items = line[1:-1].split('\t')
                geneid = int(items[0][7:])
                if (accepted is None or geneid in accepted) and geneid not in touched:
                    symbol = items[2]
                    alignment = items[3]
                    try:
                        chrm, ori, start, end = alignment.split(':')
                    except:
                        chrm, ori, start, end = '?', '?', -1,-1
                        pass
                    profile = Profile(geneid, symbol, (chrm, ori, int(start), int(end)))
                    touched.add(geneid)
                    #profiles[geneid] = profile
                else:
                    profile = None
                pass
            elif line.startswith('/'):
                if profile is not None: yield profile
                profile = None
                pass
            elif profile is not None:
                items = line.strip().split('\t')
                if items[0].startswith('A'): items = items[1:]
                position = int(items[0])
                value = float(items[2])
                error = None if len(items) < 4 else float(items[3])
                profile.set_data(position, value, error)
                pass
            pass
        pass

    @staticmethod
    def get_cpginfo(chromosome, start, end, offset):
        """drawメソッドで仕様するためのCpG islandの情報を提供するメソッド
mm8のデータを元にCpG islandか否かを配列で返す．
開始点がTSSからどれくらいずれているかoffsetに設定する"""
        import tempfile
        fn_temp = tempfile.mktemp('.txt')
        command = '/Users/takaho/Projects/C++/CpGDetector/wgcpg get -i /Data/M_musculus/NCBI36/mm_ref_chr{0}.fa -s {1} -e {2} > {3}'.format(chromosome, start, end, fn_temp)
    #print(command)
        os.system(command)
        pattern = open(fn_temp).readline().strip()
    #print(pattern)
        os.unlink(fn_temp)
        current_pat = '-'
        info = []
        start = None
        for i in range(len(pattern)):
            c = pattern[i]
            if c != current_pat:
                if c == '+':
                    start = i
                else:
                    end = i
                    if start != None: info.append((start + offset, end + offset))
                    start = None
                    pass
                current_pat = c
                pass
            pass
        if start != None: info.append((start + offset, len(pattern) + offset))
        return info
    pass


if __name__ == '__main__old__':
    accepted = []
    genes = 'Cdx2', 'Gata6', 'Hoxa5', 'Olig2', 'Adcy7', 'Cxcl14', 'Guca1a', 'Nnat'
    genes = 'Zic1', 'Hoxd11', 'Gata6', 'Cdx2'
    genes = 'Cdx2', 'Gata6', 'Hoxa5', 'Olig2'
    # similar to no antibody Cbx3, Etf1, Lias
    genes = 'Gpi1', 'Dedx31', 'Dtl', 'Eif1a', 'Eif3a','Cog4', 'Atf5', 'Hdlbp', 'Lias', 'Zfp62', 'Galk2', 'Tbp'
    #genes = 'Eif3a', 'Timm8b', 'Urb2', 'Sep15', 'Fgfr1op2', 'Fgfr1op2', 'Myl4', 'Kin', 'Rbpj', 'Parp2', 'Gpi1'
    genes = 'Zfand6', 'Orc6l', 'Ddx31', 'Rpl35a', 'Nudt6','Lig4', 'Invs','Urb2', 'Dyrk1c', 'Ints7', 'Grpel2', 'Dtl', 'Vps35', 'Rad17', 'Xrcc4', 'Pcnp', 'Cd68', 'Fkrp', 'Uchl4'
    genes = "G530011O06Rik","Fam173b","Tmem167","Nnat","8430406I07Rik","Secisbp2","Zranb1","Arhgap5","Sfrs3","Yy1","Kin","Ddx31","Grpel2","Tcte2","Matr3","Lcor","Ash1l","Dnajb4","Cd68","Klhl9","Dyrk1a","Dyrk1c","BC062115",
    genes = "Zfp180","Ssr4","Mastl","Lias","Parg","Cirh1a","Phf5a","Xrcc4","Zwilch","Ccdc127","Sept2","Lmln","Ogfod1","A630007B06Rik","Hs2st1","Snapc5","Baz2a","Plagl1","2310003L22Rik","Atp5c1","Eif1a","BC065397","Eif1ay","Zfp509","Rpl4","Cbx3","Gart","Lig4","Atp6v1d"
    genes ="Parp2","Nudt6","Rmi1","Lrrc40","Galk2","Anapc10","Paics","Mpdu1","Rbm12","Sfrs11","Pcnp","D11Wsu47e","Serpini1","Fam135a","Rhbdd3","Sep15","Mrps18b","Cog4","B230317F23Rik","Invs","Myl4","Myo18a","Mrpl47","Atp6v1a","Top3a","Xpnpep3","Ints7","Ndufs1","Polr1a"
    genes = "Ftl1","Polg2","Derl2","Rad17","Gla","Ifna13","Ptpn4","Nudcd1","Hnrnpa2b1","Hspe1","Atf5","C330027C09Rik","Hdlbp","4932432K03Rik","Eif3a","Dap3","Gpi1","Trappc4","Tbp","2410091C18Rik","Dtl"
    genes = "Gtf3c4","Srp54a","Zfand6","Nkiras1","1700040I03Rik","4930579K19Rik","Fgfr1op2","Ccdc47","2610019A05Rik","Brctd1","Zfp62","Mga","BC031181","AI316807","Fam114a2","Vps35","Rbpj","Timm8b","Orc6l"
    genes = "Zdhhc6","Mis12","Terf2ip","Nr2c2","Urb2","1500016O10Rik","Sumf2","Safb2","Tmem68","Fkrp","Iscu","Mak16","Fam162a","Rpl35a","Trove2","1110059G10Rik","Polr2b","Fam179b","1110059E24Rik","Rbm4","Uchl4","Mrpl48","Refbp2","4930583K01Rik","BC057893","Etf1"

    # available gard
    genes = 'Cdx2', 'Gata6', 'Hoxa5', 'Olig2', 'Cxcl14', 'Adcy7', 'Gart', 'Tubb1'
    genes = 'Pax5', 'Ebf1', 'Irf4', 'Egr2', 'Egr3'
    genes = 'Sox2', 'Satb2', 'Olig3', 'Tbx3', 'Olig1', 'Gata4', 'Epas1', 'Foxd3', 'Sall4', 'Id1', 'Pax5', 'Ebf1', 'Gata6', 'Cdx2'
    #genes = 'Plxnd1', 'Lhx5', 'Hoxa3', 'Hoxd13', 'Fzd1', 'Neurod2', 'Irx1', 'Vgll2', 'Cebpa', 'Pitx1', 'Shh', 'Wnt5a', 'Tbx1', 'Bmi1', 'Hoxc5', 'Hoxc6', 'Pax5', 'Satb1', 'Eya4'

    genes = 'Zic1', 'Hoxd11', 'Gata6', 'Cdx2', 'Pax5', 'Cdx2', 'Gata6', 'Hoxa5', 'Olig2', 'Adcy7', 'Cxcl14', 'Guca1a', 'Nnat'


    genes = 'Hoxd11', 'Pax3'

    genes = "Usp44","623242","Rnf220","Rnf128","Kcnc4","Zfyve28","AW146020","Limch1","St5","EG626175","Rfc3","Adra2b","ENSMUSG00000073403","Spred1","Klf15","EG328479","Bbx","Prdm5","Cdh22","Adamts15","A130092J06Rik","Rgl3","Rprml","C230071H18Rik","Axin2","EG432907","Csrnp3","1700001O22Rik","Foxp4","Slc30a2","Ccdc3","4933403G14Rik","Gpr120","She","OTTMUSG00000010975","Phpt1","Aard","Pcdh1","Rimbp3","Npm1","Msrb2","Vwa2","Atp1b2","Stk31","Atp7b","Oprd1","Bmp7","Calca","Cckbr","EG667815","EG667820","Cdh4","Cdk5r2","Cdk6","Col8a2","Chat","Rassf10","Il28ra","Pdx1","EG668855","OTTMUSG00000016072","EG668948","Epas1","Tmem132d","Gja1","EG244189","H2-Ab1","H2-D4","H2-Q1","H2-Q2","Mynn","Car15","Hes2","6430514L14Rik","Trp53i11","Frat2","Irf4","Kcna5","Kcnh3","Mto1","Kit","Krt18","Aff3","Arhgef15","Ldb2","Prokr2","Ltk","Ly75","Msc","Mapt","Neurod1","Nefm","Ngfr","Rtn4rl1","Cnnm1","Npas1","Olig1","Slc25a2","Pcsk2","Plcb2","Prph","Ptpru","Fgf12","AI593442","Rfx2","Agap2","Sfrp2","Sema3f","St8sia3","Stac2","544928","Spnb1","Npr1","Sry","Tulp1","Tyro3","Ucn","D1Pas1","Zbtb7b","Extl1","Gm1568","Clic6","Foxo3","Sez6l","Nespas","Jph3","Zc3h8","Slc18a3","6332401O19Rik","Fxyd7","EG666438","Shmt2","Ppp1r9b","Nr5a2","Psg16","Hunk","Znrf4","Gm996","Ophn1","S1pr5","Cabp7","Abcg4","Olfr441","Ocln","Garnl4","Pcsk1n","Ttc30a2","Cables1","Ap3b1","Alg11","Gabbr1","Ppp1r16b","Thpo"

    #genes = 'Tmed4','Ddx56','Armc5','Nr4a1','Tbc1d14','Olig2','Mthfd1l','Rpl37','Vasp','Stub1','Rhbdl1','Crlf1','Tmem59l','AA409316','Cd24a','Gjb3','Foxb2','Ccnd2','Tbx2','Plcb3','Ppp1r14b','Brwd1','2610017I09Rik','Pou3f3','Zbtb7b','Onecut1','Tmem204','Ndrg2','Tppp2','Syt6','BC051227','Sqle','Vsig2','Suz12','Grin2c','Nxph3','Pcdhgc4','Grin3b','Zcchc14','Syt3'

    #genes = 'AA409316','Rhbdl1','D330022A01Rik','Brwd1','Ndrg2','Tppp2','Hdac10','Grin3b','Tcfap2a','Smc6','Gen1','Rs1','Grin2c','Hspa8','EG606511','Rgl2','Cnot6','Utp15','Hoxb5','Tusc4','C1ql1','EG665291','Grin2a','432576','Kctd11','Tlcd1','Nek8','Cacna1h','Map4k5','Pkm2','Trim27','Prss34','Rela','Rpain','Pdia2','Arhgdig','Rgs11','Ercc8','Mdc1','Tubb5'
    genes = 'BC051227','Armc5','Pou3f3','Zcchc14','Lag3','A230083G16Rik','Ptms','Adamts15','Stub1','Rhbdl1','Ccnd2','LOC665268','Islr2','Lrrc56','Hras1','Padi2','1700029I15Rik','Ccdc106','Pex16','Gyltl1b','Plcb3','Ppp1r14b','Josd2','0610012D14Rik','Idh3b','Ebf4','Nkx6-2','Adrm1','Snrpb2','Grin3b','Pten','EG666410','Gna11','Hdac10','Zfp64','Aqp6','Pcdhgc5','D330022A01Rik','Shkbp1','Dus3l'
    genes = 'BC051227','Pou3f3','Armc5','Adamts15','Zcchc14','Stub1','Rhbdl1','LOC665268','Islr2','Lag3','A230083G16Rik','Ptms','Ccnd2','Lrrc56','Hras1','Josd2','0610012D14Rik','Padi2','Pten','Plcb3','Ppp1r14b','Idh3b','Ebf4','Snrpb2','Ccdc106','Nkx6-2','1700029I15Rik','Zfp64','Pex16','Gyltl1b','Aqp6','EG666410','Gna11','Adrm1','Dus3l','Pcdhgc5','D330022A01Rik','Grin3b','Shkbp1','Ubl4b'

    genes = 'Hoxa9', 'Hoxd11', 'Pax3', 'Tbx3', 'Zic1', 'Actb'

    genes = ('Ddx4', 'Rhox6', 'Mael', 'Mlf1','Dazl','Tex11', 'Irx5',
             'Hoxb4','Gata6','Eomes','Sim1','Olig2','Sox21')

    genes = ('Dazl', 'Mov10l1', 'Ddx4', 'Lef1', 'Bmp6', 'Fbxo46', '2900010M23Rik', 'Hdac2', 'Insm1', 'Cdx2', 'Ccnd2', 'Dll4', 'Igf2', 'Slc35d3', 'Pcdh7')
#    genes = ('Dazl', 'Bmp6')

    if 0:
        genes = {}
        import re, os, sys
        import mousegene
        for gene in mousegene.MouseGene.get_all():
            m = re.search('Hox(a|b|c|d)(\\d+)$', gene.symbol)
            if m:
                genes[gene.symbol] = m.group(1) + ('00000{0}'.format(100 - int(m.group(2))))[-4:]
                #print(gene.symbol, genes[gene.symbol])
                pass
            pass
        genes = sorted(genes.keys(), key=lambda g_: genes[g_])
        pass

    genes = ('Six2', 'Vsx2', 'Hist3h2a', 'Hist3h2bb', 'Rnf187', 'Tcf21', 'Nkx6-2', 'Dmbx1', 'Slc32a1', 'Isl2', 'Ccdc108', 'Ralgps2', 'Barhl1', 'Gad1', 'Wt1', 'Shisa3', 'Gsx2', 'Cux2', 'Evx1', 'Foxa1', 'Barx1', 'Ankrd34b', 'EGG667363', 'Mafa', 'Adcy5')

    genes = ('Evx', 'Hoxd13', 'Gata4', 'Pax3', 'Hoxa10', 'Otx2', 'Hand2', "Pax6",
             'Kank1', 'Ccr10', 'Chd5', 'Arpc5', 'Gfi1b')
    genes = 'Hoxa9', 'Pax9', 'Tbx3', 'Jak2', 'Ptx3', 'Gsn', 'Itga2', 'Tuba1a'


    genes = ('Fndc3c1','Ipw''Oasl1',
             'Sox17','Bmi1','Nkx2-6','Pde10a','Tbr1','Tbr4')

    import mapview
    seq_gene = mapview.get_seq_gene()
    s2g = {}
    for gene in seq_gene.values():
        s2g[gene.symbol] = gene
        pass

    for symbol in genes:
        gene = s2g.get(symbol, None)
        #gene = mousegene.get_gene_by_name(symbol)
        if gene is not None:
            accepted.append(gene.geneid)
        else:
            sys.stderr.write('%s does not have gene\n' % symbol)
            pass
        pass
    #accepted = set(accepted)
    import reportlab.pdfgen.canvas
    cpgpatterns = {}
    colors = 'green', 'navy', 'crimson', 'olive', 'cadetblue', 'brown'
    colors = 'navy', 'red', 'olive', 'orange', 'green', 'crimson', 'blue', 'khaki', 'cadetblue', 'pink',
    colors = 'ReportLabGreen', 'navy', 'brown'
    left, top, right, bottom = 100, 750, 600, 50
    gx, gy, gw, gh = left, top, 120, 100

    filenames = []
    args = sys.argv[1:]
    args.reverse()
    use_cpg = False
    filename_output = 'profiles.pdf'
    verbose = False
    tss_left = -6000
    tss_right = 4000
    multiplier = 1.0
    while len(args) > 0:
        arg = args.pop()
        if os.path.exists(arg) and arg.endswith('.txt') and arg.startswith('-') is False:
            filenames.append(arg)
        elif arg == '-cpg':
            use_cpg = True
        elif arg == '-o':
            filename_output = args.pop()
        elif arg == '-verbose':
            verbose = True
        elif arg == '-ymax':
            Profile.MAXIMUM_VALUE = 2 ** int(args.pop())
        elif arg == '-left':
            tss_left = int(args.pop())
        elif arg == '-right':
            tss_right = int(args.pop())
        elif arg == '-m':
            multiplier = float(args.pop())
            pass
        pass
    #filenames.reverse()

    pagesize = (gw + 80) * len(genes) + 200, (len(filenames) * (gh + 40) + 200)
    top = pagesize[1] - gh - 60
    cnv = reportlab.pdfgen.canvas.Canvas(filename_output, pagesize=pagesize)
    y = top
    #Profile.MAXIMUM_VALUE=2 ** 4
    Profile.set_range(tss_left, tss_right)

    for i, filename in enumerate(filenames):
        if verbose: sys.stderr.write('{0}\t{1}\n'.format(i + 1, filename))
        bar_color = colors[i % len(colors)]
        name = filename.split('/')[-1].split('.')[0]
        x = left
        cnv.drawString(0, y, name)
        for prof in Profile.read_profile(filename, set(accepted)):
            for index, geneid in enumerate(accepted):
                if geneid == prof.geneid:
                    x = left + (gw + 80) * index
                    break
                pass

            cpginfo = None
            if use_cpg:
                cpginfo = []
                import os, sys
                import cpg
                import tempfile
                if prof.geneid in cpgpatterns:
                    pattern = cpgpatterns[prof.geneid]
                else:
                    filename = ' ~/Projects/C++/CpGDetector/wgcpg get -verbose -i /Data/M_musculus/NCBI36.1/mm_ref_chr5.fa -s 147615597 -e 147619597'
                    al = prof.alignment.split(':')
                    if al[1] == '-':
                        left_pos = int(al[3]) - 4000
                        right_pos = left_pos + 10000
                    else:
                        left_pos = int(al[2]) - 6000
                        right_pos = left_pos + 10000
                        pass
                    fn_temp = tempfile.mktemp('.txt')
                    cmd = '/Users/takaho/Projects/C++/CpGDetector/wgcpg get -verbose -i /Data/M_musculus/NCBI36.1/mm_ref_{0}.fa -s {1} -e {2} > {3}'.format(al[0], left_pos, right_pos, fn_temp)
                    sys.stderr.write(cmd + '\n')
                    os.system(cmd)
                    pattern = open(fn_temp).readlines()[-1]
                    os.unlink(fn_temp)
                    if al[1] == '-':
                        p_ = ''
                        for i in range(len(pattern)):
                            p_ += pattern[-1 - i]
                            pass
                        pattern = p_
                        pass
                    cpgpatterns[prof.geneid] = pattern
                    pass
                start = end = -1
                for i in range(len(pattern)):
                    p = pattern[i]
                    if p == '-':
                        if start >= 0:
                            end = i
                            cpginfo.append((start - 6000, end - 6000))
                            start = -1
                            pass
                        pass
                    elif p == '+':
                        if start < 0: start = i
                        pass
                    pass
                if start >= 0: cpginfo.append((start - 6000, end - 6000))
                pass
            prof.draw(cnv, (x, y, x + gw, y - gh), bar_color=bar_color, bar_width=0.5, cpginfo=cpginfo)
            #x += gw + 200
            pass
        y -= gh + 40
        pass
    cnv.save()
    pass

def retrieve_profiles(sample_name, chromosome, start, end, directory='mm9'):
    """read .dat of Takaho's chip-seq data"""
    ft = tempfile.mktemp('.txt')
    cmd = '/Users/takaho/Projects/C++/ChIPSeq/decseq -d "{0}" -n "{1}" -p {2}:{3}:{4} > {5}'.format(directory, sample_name, chromosome, start, end, ft)
    #print(cmd)
    os.system(cmd)
    data = {}
    for line in open(ft):
        items = line.strip().split('\t')
        if len(items) >= 3:
            data[int(items[0])] = int(items[1])
            pass
        pass
    os.unlink(ft)
    return data

def parzen(dist, size):
    #print(dist, size)
    x = abs(float(dist)) * 2 / size
    if x > 2.0: return 0.0
    if x > 1.0:
        return (2.0 - x) ** 3 * 0.36363867770067627703
    else:
        return (1.0 - 1.5 * x ** 2 + 0.75 * x ** 3) * 1.45455471080270510812;
    pass

def get_series(data, start=None, end=None, division=None, smoothlevel=None):
    if len(data) == 0:
        sys.stderr.write('empty data\n')
        return []
        pass
    if start is None: start = min(data.keys())
    prev = None
    d_min = None
    for s in sorted(data.keys()):
        if start is None: start = s
        if prev is not None:
            d = s - prev
            if d_min is None or d < d_min: d_min = d
            pass
        prev = s
        pass
    if division is None: division = d_min
    if end is None: end = prev
    #end = prev
    #division = d_min
    num_items = (end - start + division - 1) // division #+ smoothlevel
    #if end is not None: end = num_items
    values = [0] * num_items
    for p, val in data.items():
        n = (p - start) // division
        #print(n, num_items, p)
        values[n] = val
        pass
    if smoothlevel:
        coeff = []
        for i in range(smoothlevel):
            coeff.append(parzen(i, smoothlevel))
            pass
        #print(coeff)
        smoothened = []
        for i in range(num_items):
            j = - smoothlevel + 1
            w = 0.0
            v = 0.0
            while j < smoothlevel:
                if 0 <= (i+j) < num_items:
                    c = coeff[abs(j)]
                    w += c
                    v += c * values[i+j]
                    pass
                j += 1
                pass
            if w > 0.0:
                smoothened.append(v / w)
            else:
                smoothened.append(0.0)
                pass
            pass
        values = smoothened
        pass
    return values

def retrieve_cgi_ucsc(chromosome, start, end, filename='/Data/Genomes/M_musculus/mm9/cpgIslandExt.txt'):
    if chromosome.startswith('chr'):# is False:
        chromosome = chromosome[3:]#'chr{0}'.format(chromosome)
        pass
    location = []
    for line in open(filename):
        items = line.split('\t')
        if items[1][3:] == chromosome:
            p = int(items[2])
            q = int(items[3])
            if p <= end and q >= start:
                print(chromosome, start, end, p, q)
                location.append([p,q])
                pass
            pass
        pass
    return location


def draw_chipseq_profiles(genes, experiments, basedir='mm9', canvas=None, normalizer=None):
    bottom = 650
    box_height = 50
    ystep = 70
    if canvas is None:
        page_height = (len(experiments) * (len(genes) + 1) + 1) * ystep
        canvas = reportlab.pdfgen.canvas.Canvas('canvas.pdf', pagesize=(500,page_height) )
        bottom = page_height - ystep - 30
        pass
    print('drawing')
    ymax = 0.25#15#1.0#0.20#1.0#0.15
    ymax = 0.10#15#1.0#0.20#1.0#0.15
    margin=6000
    #ymax = 1.0
    #ymax = 0.25
    coefficient = 0.00005
    coefficient = 0.00001
    #coefficient = 0.00005

    # 91009101
    # 66935806
    # 80965679
    # 75417405
    #mouse_gene = mapview.get_seq_gene()
    #s2g = {}
    #for gene in genes:#mouse_gene.values():
    #s2g[gene.symbol] = gene
    #pass

    small_arrow = canvas.beginPath()
    small_arrow.moveTo(0,0)
    small_arrow.lineTo(0,2)
    small_arrow.lineTo(4,0)
    small_arrow.lineTo(0,-2)
    small_arrow.lineTo(0,0)
    small_arrow.close()

    colors = ['ReportLabGreen', 'olive', 'ReportLabFidBlue', 'cadetblue', 'ReportLabFidRed', 'tomato', 'gray', 'purple', 'brown', 'navy', 'gold', 'pink']

    for gene in genes:#symbol in symbols:#selected:
        print(gene.symbol)
        #gene = s2g[symbol]
        symbol = gene.symbol
        chrm = gene.chromosome
        ori = gene.orientation
        start = gene.position5
        end = gene.position3

#    for loc in data.split('\n'):
        #symbol, chrm, ori, start, end = loc.split(',')
        left = ((int(start) - margin) // 1000) * 1000
        right = ((int(end) + margin + 999) // 1000) * 1000
        #print(10000.0 / (right - left))
        canvas.setFillColor('black')
        canvas.setFont('Helvetica', 10)
        canvas.setFillColorRGB(0.02, 0.01, 0.02)
        canvas.drawString(100, bottom + box_height + 10, symbol)

        for index, sample in enumerate(experiments):
            #print(sample)
            if sample.find('/') >= 0:
                pos = sample.rfind('/')
                directory, sample = sample[0:pos], sample[pos + 1:]
            else:
                directory, sample = basedir, sample
                pass
            prof = retrieve_profiles(sample, chrm, left, right, directory)#basedir)
#            prof = retrieve_profiles(sample, chrm, left, right, basedir)
#                                      '/Users/takaho/SkyDrive/Research/mm9_hase/')
            #coefficient = normalizer[index]
            xconv = lambda pos: float(pos - left) * coefficient * 300 + 100
            #xconv = lambda pos: float(pos - left) / (right - left) * 300 + 100
            ymax_ = ymax #* normalizer[index]
            yconv = lambda val: min(ymax_, val) / (ymax_) * box_height + bottom
            if sample == 'Treg_WT_':
        #scale
                canvas.setStrokeColorRGB(.02, .04, .01)
                xs, ys = xconv(right) - 40, bottom + 65
                bar_width = xconv(1000) - xconv(0)
                canvas.line(xs, ys, xs - bar_width, ys)
                canvas.setFillColor('black')
                canvas.drawString(xs + 8, ys - 4, '{0} kb'.format(1))
                pass
            x0 = xconv(left)
            x1 = xconv(right)
            y0 = yconv(0)
            y1 = yconv(ymax)# * normalizer[index])
            div = 25
            prof = get_series(prof, start=left, end=right, division=div, smoothlevel=8)
            path = None
            for i, v in enumerate(prof):
                x = xconv(i * div + left)
                #print(i, v)
                y = yconv(v * 0.01)# * normalizer[index])
                #y = yconv(v * normalizer[index])
                if path is None:
                    path = canvas.beginPath()
                    path.moveTo(x, y0)
                    pass
                path.lineTo(x, y)
                pass
            if path is not None:
                path.lineTo(x, y0)
                path.close()
                canvas.setFillColor(colors[index % len(colors)])#'ReportLabFidBlue')
                canvas.drawPath(path, 0, 1)
                pass
            canvas.saveState()
            canvas.setStrokeColor('black')
            canvas.rect(x0, y0, x1 - x0, y1 - y0)

            canvas.rect(x0, y0, x1 - x0, y1 - y0)
            canvas.setFillColorRGB(.01, .0, .0)
            canvas.drawString(x0 - 100, y1 - 12, sample.replace('_', ' '))
            canvas.restoreState()

            bottom -= ystep
            pass

        canvas.setFillColorRGB(.2,.2,.2)
        canvas.setFont('Courier', 6)
        for p in gene.position5 - margin, gene.position3 + margin:
            if p < gene.position5:
                label = 'chr{0}:{1}'.format(gene.chromosome, p)
            else:
                label = '{0}'.format(p)
                pass
            canvas.drawString(xconv(p) - canvas.stringWidth(label)/ 2, bottom + 60 , label)
            pass
        y = bottom + 40
        first = True
        for rna in sorted(gene.rna, key=lambda r:r.end - r.start, reverse=True):
            x0 = xconv(rna.start)
            x1 = xconv(rna.end)
            if first:
                path = canvas.beginPath()
                if ori == '+':
                    x2, x3, x4 = x0, x0 + 20, x0 + 17
                else:
                    x2, x3, x4 = x1, x1 - 20, x1 - 17
                    pass
                path.moveTo(x2, y + 3)
                path.lineTo(x2, y + 8)
                path.lineTo(x3, y + 8)
                path.moveTo(x4, y + 11)
                path.lineTo(x3, y + 8)
                path.lineTo(x4, y + 5)
                canvas.setStrokeColor('olive')
                canvas.drawPath(path)
                first = False
                pass
            canvas.setStrokeColor('gray')
            canvas.line(x0, y, x1, y)
            xa = x0 + 5
            while xa + 10 < x1:
                canvas.saveState()
                canvas.setFillColor('gray')
                canvas.translate(xa, y)
                if ori == '-':
                    canvas.rotate(180)
                    pass
                canvas.drawPath(small_arrow, 0, 1)
                canvas.restoreState()
                xa += 25
                pass

            for exon in rna.exons:
                if exon[0] == 'CDS':
                    canvas.setFillColor('crimson')
                else:
                    canvas.setFillColor('cadetblue')
                    pass
                x0 = xconv(exon[1])
                x1 = xconv(exon[2])
                canvas.rect(x0, y - 3, x1 - x0, 6, 0, 1)
                pass
            y -= 8
            break
        bottom -= ystep
        pass
    return canvas

class ProfileData(object):
    def __init__(self, pos=None, value=None):
        if pos is not None and value is not None:
            self.__position = pos
            self.__value = value
            self.__count = 1
        else:
            self.__position = self.__value = self.__count = 0
    def append(self, value):
        self.__value += value
        self.__count += 1

    position = property(lambda s:s.__position)
    count = property(lambda s:s.__count)
    value = property(lambda s:s.__value)

def get_average_profile(filename, division=100, upstream=-10000, downstream=10000, genes=None, cache=None):
    import sqlite3, pickle, md5, base64, zlib

    if cache is None:
        cache = os.path.expanduser('~/.chipprofile.db')

    if genes is None:
        genekey = ''
    else:
        import md5
        m = md5.md5()
        for gene in sorted(genes):
            m.update(gene)
        genekey = m.hexdigest()
    with sqlite3.Connection(cache) as cnx:
        data = None
        cur = cnx.cursor()
        cur.execute('select name from sqlite_master where name=?', ('tssprofs', ))
        r = cur.fetchone()
        if r is None:
            #sys.stderr.write('create table\n')
            cur.execute('create table tssprofs (filename not null, upstream int4 not null, downstream int4 not null, division int4, genes, data blob)')
            cnx.commit()
        else:
            #sys.stderr.write('load data\n')
            cur.execute('select data from tssprofs where filename=? and upstream <= ? and downstream >= ? and genes=? and division=?', (filename, upstream, downstream, genekey, division))
            r = cur.fetchone()
            #print(r)
#            print(genekey, filename, upstream, downstream)
            if r is not None:
                data = pickle.loads(zlib.decompress(base64.b64decode(r[0])))
        if data is None:
            touched = set()
            data = {}
            line_num = 0
            with open(filename) as fi:
                line_num += 1
                name = None
                for line in fi:
                    if line.startswith('>'):
                        items = line.strip().split('\t')
                        name = items[2]
                        if name in touched or (genes is not None and name not in touched):
                            name = None
                        else:
                            touched.add(name)
                            sys.stderr.write(' {} {}              \r'.format(len(touched), name))
                    elif line.startswith('/'):
                        name = None
                    elif name is not None:
                        items = line.split('\t')
                        try:
                            pos = int(items[0])
                            value = float(items[2])
                        except:
                            print(items)
                            print(line_num)
                            print(filename)
                            raise
                        if pos not in data:
                            data[pos] = ProfileData(pos, value)
                        else:
                            data[pos].append(value)
            cur.execute('insert into tssprofs values(?, ?, ?, ?, ?, ?)', (filename, upstream, downstream, division, genekey, base64.b64encode(zlib.compress(pickle.dumps(data)))))
            sys.stderr.write('committed {}\n'.format(filename))
            cnx.commit()
    N = (downstream - upstream + division - 1) // division
    accum = [0] * N
    counts = [0] * N
    for pos in sorted(data.keys()):
        section = int(math.floor(float(pos - upstream) / division))
        if 0 <= section < N:
            val = data[pos]
            accum[section] += val.value
            counts[section] += val.count
    values = []
    for i in range(N):
        values.append(accum[i] / counts[i] if counts[i] > 0 else 0)
    #print(values)
    return values

def draw_total_profiles():
    import argparse, collections
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', default='out.pdf', help="output filename")
    parser.add_argument('-i', nargs='+', required=True, help="show this help")
    #parser.add_argument('--direcotry', help='src directory', default='mm9')
    #parser.add_argument('-y', default=None)
    parser.add_argument('-u', type=int, default=-10000, help="upstream distance")
    parser.add_argument('-d', type=int, default=10000, help="downstream distance")
    parser.add_argument('-s', type=int, default=100, help="step")
    parser.add_argument('-w', type=int, default=100, help="window size")
    parser.add_argument('--smooth', type=int, default=0, help="smoohing level (default:0)")
    parser.add_argument('-g', default=None, help="gene list")
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('-y', default=None, help="Maximum value of Y axis")
    args = parser.parse_args()

    if args.y:
        ymin, ymax = [float(x_) for x_ in args.y.split(',')]
    else:
        ymin = 0
        ymax = None
    verbose = args.verbose
    data = collections.OrderedDict()

    genes = None
    upstream = args.u
    downstream = args.d
    window_size = args.w

    if args.g:
        genes = set()
        with open(args.g) as fi:
            for line in fi:
                items = line.strip().split('\t')
                genes.add(items[0])

    for filename in args.i:
        name = os.path.basename(filename)
        rpos = name.rfind('.')
        if rpos > 0: name = name[0:rpos]
        if verbose:
            sys.stderr.write('{}\n'.format(filename))
        data[name] = get_average_profile(filename, window_size, upstream, downstream, genes)

    #gr = reportlab.pdfgen.canvas.Canvas(args.o)
    pos = upstream

    line = '#pos\t' + '\t'.join(data.keys())
    if ymax is None:
        ymax = 0
        for val in data.values():
            ymax = max(ymax, max(val))

    import tkgraph
    gr = tkgraph.LineGraph(args.o, tkgraph.Axis(upstream, downstream, 50, 200, major_tick=max(1000, -upstream)),
                       tkgraph.Axis(ymin, ymax, 500, 650, major_tick=ymax - ymin))

    ycnv = lambda y:min(y, ymax) * height + bottom
    colors = ('ReportLabFidRed', 'ReportLabGreen', 'ReportLabFidBlue', 'cadetblue', 'violet', 'pink', 'gold', 'olive', 'cyan', 'magenta', 'crimson', 'navy', 'gray', 'purple', 'brown', )

    index = 0
    for name, beddata in marks:
        for region in gtf.Section.find(beddata,
    for name, values in data.items():
        color = colors[index % len(colors)]
        gr.begin_item(name)
        path = None
        for i, v in enumerate(values):
            x = gr.x(i * window_size + upstream)
            y = gr.y(v)
            if path is None:
                path = gr.beginPath()
                path.moveTo(x, y)
            else:
                path.lineTo(x, y)
        gr.setStrokeColor(color)
        gr.drawPath(path)
        ly = gr.top - index * 10
        gr.line(gr.right + 100, ly, gr.right + 120, ly)
        gr.drawString(gr.right + 130, ly, name)
        index += 1
    gr.save()

    print(line)
    index = 0
    while pos <= downstream:
        line = '{}'.format(pos)
        for val in data.values():
            if index < len(val):
                line += '\t{:.3f}'.format(val[index])
            else:
                line += '\t.'
        print(line)
        pos += window_size
        index += 1

    #gr.save()


if __name__ == '__main__':

    draw_total_profiles()
    exit()

    import pickle
    print('loading')
#    genes = pickle.load(open('.cache'))
    genes = ('Foxo1', 'Ebf1', 'Fosb', 'Pax5', 'Ikzf3', )
    experiments = ('S02711', 'S02712', 'S02713')
    basedir = 'tss'
    draw_chipseq_profiles(genes.values(), experiments, basedir=basedir).save()
    exit()

    experiments = ['Ring1b-WTMEF_', 'Ring1b-WTMEF1_', 'Ring1b-L307RMEF_', 'Ring1b-L307RMEF1_', 'H3K27me3-MEFWT_', 'H3K27me3-MEFKO_', '/Data/ChIPSeq/Mendoh/RING1B/mm9/Ring1B_R1ABdKO_']
    print('try to draw')

    draw_chipseq_profiles(genes.values(), experiments, basedir='/Volumes/ThunderDrive/Isono/ChIP-seq Isono/nd/').save()
    pass
