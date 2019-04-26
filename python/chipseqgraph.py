#coding:utf-8

import os, sys, re, subprocess, math, argparse, sqlite3, gzip, pickle
import reportlab.pdfgen.canvas
from genomicitem import GenomicItem as GI

__path = '/Users/takaho/Projects/C++/ChIPSeq/'
bedcache = {}
def load_peaks(filename, chrom, start, stop):
    if os.path.exists(filename) is False: return []
    margin = 5000
    peaks = []
    with open(filename) as fi:
        for line in fi:
            items = line.strip().split('\t')
            if line.startswith('#') or items[1].isdigit() is False: continue
            p5 = int(items[1])
            p3 = int(items[2])
            if items[0] == chrom and p5 < stop + margin and p3 > start - margin:
                peaks.append([p5, p3])
    return sorted(peaks, cmp=GI.compare_position)

def load_bed_sections(filename, chromosome, start, stop):
    global bedcache
    #print(filename)
    if filename in bedcache:
        peaks = bedcache[filename]
        for peak in GI.find_in_range(peaks, chromosome, start, stop):
            yield peak
    else:
        peaks = []
        #print(chromosome, start, stop)
        with open(filename) as fi:
            for line in fi:
                if line.startswith('#'): continue
                items = line.strip().split('\t')
    #            print(items)
                if len(items) >= 3 and items[1].isdigit() and items[2].isdigit():
                    if re.match('chr[0-9A-Z]+\\s', line):
                        ch = items[0]
                        p1 = int(items[1])
                        p2 = int(items[2])
                        peak = GI(ch, '+', p1, p2)

                        if chrom == ch and p2 >= start and p1 <= stop:
                            #print(line)
                            yield peak
                        peaks.append(peak)
        bedcache[filename] = sorted(peaks, cmp=GI.compare_position)#peaks



def get_profile(basedir, name, chrom, start, stop, coeff=1.0, division=25):
    if chrom.startswith('chr') is False: chrom = 'chr' + chrom
    cmd = os.path.join(__path, 'decseq'), '-d', basedir, '-n', name, '-p', '{}:{}:{}'.format(chrom, start, stop), '-w', '{}'.format(division)
    #print(' '.join(cmd))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    data = {}
    for line in proc.stdout:
        items = line.strip().split('\t')
        if line.startswith('>') or line.startswith('/'): continue
        items = line.strip().split('\t')
        data[int(items[0])] = int(items[1]) * coeff
    return data

_chipchipdata = {}
def get_chipchip_profile(filename, chromosome, start, stop):
    global _chipchipdata
    if filename in _chipchipdata:
        chipchip = _chipchipdata[filename]
        data = chipchip[chromosome]
    else:
        chipchip = {}
        #print(filename)
        if filename.endswith('.wig'):
            with open(filename) as fi:
                fi.readline()
                chrom = None
                fixed = False
                for line in fi:
                    #sys.stdout.write('{} {}'.format(len(chipchip), line))
                    if line.startswith('variable'):
                        m = re.search('chrom=(\\w+)', line)
                        chrom = m.group(1)
                        if chrom.startswith('chr') is False: chrom = 'chr' + chrom
                        m = re.search('span=(\\d+)', line)
                        span = int(m.group(1))
                        chipchip[chrom] = []
                        pos = 0
                        step = 0
                        fixed = False
                        print('variable ' + chrom)
                    elif line.startswith('fixedStep'):
                        m = re.search('chrom=(\\w+)', line)
                        chrom = m.group(1)
                        if chrom.startswith('chr') is False: chrom = 'chr' + chrom
                        m = re.search('step=(\\d+)', line)
                        step = int(m.group(1))
                        span = step
                        chipchip[chrom] = []
                        m = re.search('start=(\\d+)', line)
                        if m: pos = int(m.group(1))
                        fixed = True
                        #span = int(m.group(1))
                        #chipchip[chrom] = []
                        print('fixed ' + chrom)
                    elif chrom is not None:
                        if fixed:
                            chipchip[chrom].append([None, chrom, pos, pos + step, -float(line) * math.log(2)])
                            pos += step
                        else:
                            items = line.strip().split('\t')
                            pos = int(items[0])
                            value = float(items[1])
                        #print(start, span, value)
                            chipchip[chrom].append([None, chrom, pos, pos + span, -value * math.log(2)])
            _chipchipdata[filename] = chipchip
        else:
            with open(filename) as fi:
                # ? chromosome start stop value(log10)
                for line in fi:
                    items = line.strip().split('\t')
                    #print(filename, items)
                    chrom = items[1]
                    if chrom not in chipchip:
                        chipchip[chrom] = []
                    chipchip[chrom].append([items[0], items[1], int(items[2]), int(items[3]), float(items[4]), float(items[5]), float(items[6])])
        _chipchipdata[filename] = chipchip
    data = chipchip.get(chromosome, [])
    # for c in chipchip.keys():
    #     print('{}\t{}'.format(c, len(chipchip[c])))
    #     #print(chipchip[c][0:3])
    #     #chipchip[c] = sorted(chipchip[c], key=lambda x_: x_[1])
    #     #print('{}\t{}'.format(c, len(chipchip[c])))
    #     #print(chipchip[c][0:3])
    left = 0
    right = len(data)
    index = None
    while left < right:
        center = (left + right) // 2
        datum = data[center]
        s, e = datum[2], datum[3]
        #print(center, len(data), s, e, start, stop)
        if s < stop and e >= start:
            index = center
            break
        if s >= stop:
            right = center
        elif e < start:
            left = center + 1
    #
    if index is None: return {}
    i = index
    left_limit = 0
    right_limit = len(data)
    while i >= 0:
        datum = data[i]
        s, e = datum[2], datum[3]
        if e < start:
            left_limit = i
            break
        i -= 1
    i = index + 1
    while i < len(data):
        datum = data[i]
        s, e = datum[2], datum[3]
        if s > stop:
            right_limit = i - 1
            break
        i += 1
    prof = {}
    #print(left_limit, right_limit, left, right)
    for datum in data[left_limit:right_limit]:
        pos = (datum[2] + datum[3]) // 2
        #print(pos, datum[4])
        prof[pos] = -float(datum[4]) / math.log10(2)
    #print(prof)
    return prof

def load_genes(filename_gtf, chromosome, start, stop, dbfile=None):
    import gtf
    return gtf.get_genes(filename_gtf, chromosome, start, stop, dbfile)
    # table_name = 'genes'
    # if dbfile is None:
    #     dbfile = os.path.join(os.path.dirname(filename_gtf), '.{}.db'.format(os.path.basename(filename_gtf)))
    #     if os.path.exists(dbfile) and os.path.getsize(dbfile) > 10000:
    #         cnx = sqlite3.Connection(dbfile)
    #         cur = cnx.cursor()
    #         cur.execute('select * from sqlite_master where name="{}"'.format(table_name))
    #         r = cur.fetchone()
    #         if r is not None:
    #             cur.execute('select count(*) from {}'.format(table_name))
    #             r = cur.fetchone()
    #             if r is not None and r[0] > 10000:
    #                 cur.execute('select id, name, chromosome, orientation, start, stop, exons, coding in {} where chromosome="{}" and start <= {} and stop >= {} and coding=1'.format(table_name, chromosome, stop, start))
    #                 for r in cur.fetchall():
    #                     gene = Gene(r[0], r[1], r[2], r[3])
    #                     for exon in r[6].split(','):
    #                         s, e = [int(x_) for x_ in exon.split('-')]
    #                         gene.add_exon(s, e)
    #                     genes.append(gene)
    #                 cur.close()
    #                 cnx.close()
    #                 return genes
    # margin = 1000000
    # #allgenes = []
    # genes = {}
    # accepted = {}
    # coding = set()
    # with open(filename_gtf) as fi:
    #     for line in fi:
    #         items = line.strip().split('\t')
    #         p5 = int(items[3])
    #         p3 = int(items[4])
    #         #print(items[0], chromosome, p5, p3, start, stop, p5 - margin - stop)
    #         #if p3 + margin < start or p5 - margin > stop: continue
    #         info = {}
    #         for elem in items[8].split(';'):
    #             m = re.match('(\\w+)\\s+"([^"]+)', elem.strip())
    #             if m:
    #                 info[m.group(1)] = m.group(2)
    #         if 'transcript_id' not in info or 'gene_name' not in info: continue
    #         trid = info.get('transcript_id', None)
    #         if trid not in genes:
    #             gene = Gene(trid, info['gene_name'], items[0], items[6])
    #             #print(gene.name, len(allgenes))
    #             #allgenes.append(gene)
    #             #if items[0] == chromosome:
    #             genes[trid] = gene
    #         else:
    #             gene = genes[trid]
    #         if items[2] != 'exon':
    #             if items[2] == 'start_codon' or items[2] == 'stop_codon':
    #                 if trid is not None:
    #                     if trid not in accepted: accepted[trid] = 0
    #                     accepted[trid] |= 1 << (items[2] == 'stop_codon')
    #         else:
    #             gene.add_exon(p5, p3)
    # inside = []
    # for trid, flag in accepted.items():
    #     #print(trid, flag)
    #     g = genes[trid]
    #     if flag == 3 and g.chromosome == chromosome and g.start <= stop and g.stop >= start:
    #         inside.append(genes[trid])
    # cnx = sqlite3.Connection(dbfile)
    # cur = cnx.cursor()
    # cur.execute('select * from sqlite_master where name="{}"'.format(table_name))
    # r = cur.fetchone()
    # if r is None:
    #     cur.execute('create table {} (id not null primary key, name, chromosome, orientation, start int4, stop int4, exons, coding int4)'.format(table_name))
    #     cur.execute('create index index_{0} on {0} (chromosome, start, stop)'.format(table_name))
    #     pass
    # for gene in genes.values():
    #     coding = int(accepted.get(gene.transcript_id, 0) == 3)
    #     estr = ''
    #     for p1, p2 in gene.exons:
    #         if estr != '': estr += ','
    #         estr += '{}-{}'.format(p1, p2)
    #     cur.execute('insert into {} values("{}", "{}", "{}", "{}", {}, {}, "{}", {})'.format(table_name, gene.transcript_id, gene.name, gene.chromosome, gene.orientation, gene.start, gene.stop, estr, coding))
    # cur.close()
    # cnx.commit()
    # cnx.close()
    # return sorted(inside, key=lambda g:g.start)

def parzen(dist, size):
    #print(dist, size)
    x = abs(float(dist)) * 2 / size
    if x > 2.0: return 0.0
    if x > 1.0:
        return (2.0 - x) ** 3 * 0.36363867770067627703
    else:
        return (1.0 - 1.5 * x ** 2 + 0.75 * x ** 3) * 1.45455471080270510812;
    pass

def smoothen(values, smoothlevel):
    coeff = []
    for i in range(smoothlevel):
        coeff.append(parzen(i, smoothlevel))
        pass
    smoothened = []
    num_items = len(values)
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
    return smoothened

def convert_dict_to_profile(prof, start, stop, division, totalcounts=None):
    #print(totalcounts)
    smin = 0#min(counts.keys())
    smax = (stop - start) // division + 1
    values = [0] * (smax - smin)# + 1)
    for section, count in prof.items():
        sec = (section - start) // division
        if 0 <= sec < len(values):
            values[sec] += count
    if totalcounts is not None:
        coeff = 1e9 / division / totalcounts
    else:
        coeff = 1e3 / division
#    coeff *= 25.0 / division
    return [v_ * coeff for v_ in values]

def convert_dict_to_chipchip_profile(prof, start, stop, division):
    smin = 0#min(counts.keys())
    smax = (stop - start) // division + 1
    values = [0] * (smax - smin + 1)
    counts = [0] * (smax - smin + 1)
    for position, value in prof.items():
        sec = (position - start) // division
        if 0 <= sec < len(values):
            values[sec] += value
            counts[sec] += 1
    return [((values[i] / counts[i]) if counts[i] > 0 else 0.0) for i in range(len(values))]

def get_step(minval, maxval, minsteps=2, maxsteps=10):
    f = 0
    r = [1.0,2.0,5.0]
    dif = maxval - minval
    index = 0
    while 1:
        t = r[index] * (10 ** f)
        n = math.floor(dif / t)
        if minsteps <= n <= maxsteps:
            return t
        if n > maxsteps:
            index += 1
            if index >= len(r):
                index = 0
                f += 1
        if n < minsteps:
            index -= 1
            if index < 0:
                index = len(r) - 1
                f -= 1

#
def count_total_sequences(srcdir, name):
    cachefile = os.path.join(srcdir, '.' + name + '.total.cache')
    if os.path.exists(cachefile) and os.path.getsize(cachefile) > 0:
        with gzip.open(cachefile) as fi:
            data = pickle.load(fi)
            return data
    #
    cmd = ['/Users/takaho/Projects/C++/ChIPSeq/countall']
    for fn in os.listdir(srcdir):
        if fn.endswith('.dat') and fn.startswith(name + 'chr') > 0 and fn.find('_random') < 0 and fn.find('chrUn') < 0:
            cmd.append(os.path.join(srcdir, fn))
    #
    if len(cmd) < 2: return None
    sys.stderr.write(' '.join(cmd) + '\n')
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    data = {}
    for line in proc.stdout:
        items = line.strip().split('\t')
        m = re.search('(chr[0-9XYMT]+)', items[0])
        if m:
            data[m.group(1)] = int(items[1])
    #
    #print(data)
    with gzip.open(cachefile, 'wb') as fo:
        pickle.dump(data, fo)
    return data

def draw_barcode_graphs(args):
    total_counts = {}
    for name in args.n:
        total_counts[name] = count_total_sequences(args.d, name)
        #print(count_total_sequences(args.d, name))
        #print(name, sum(total_counts[name].values()))
    geo = [float(x_) for x_ in args.g.split(',')]
    width = geo[0]
    height = geo[1]
    cnv = reportlab.pdfgen.canvas.Canvas(args.o)
    division = args.w
    left = 50
    top = 700
    vmin = 0
    vmax = 10
    x = left
    spacer = 15
    for location in args.l:
        chrom, pos = location.split(':')
        start, stop = [int(x_) for x_ in pos.split('-')]
        xcnv = lambda index: float(index * division) / (stop - start) * width + x
        y = top
        for name in args.n:
            tc = sum(total_counts[name].values())
            prof = get_profile(args.d, name, chrom, start, stop)
            values = convert_dict_to_profile(prof, start, stop, division, tc)
            if args.s > 0:
                values = smoothen(values, args.s)
            #
            for i, v in enumerate(values):
                x0 = xcnv(i)
                x1 = xcnv(i + 1)
                degree = min(1.0, max(0, (v - vmin) / (vmax - vmin)))
                if name.startswith('Liver'):
                    cnv.setFillColorRGB(1.0 - degree, 0.8 + (1.0 - degree) * 0.2, 1.0 - degree)
                elif name.startswith('Heart'):
                    cnv.setFillColorRGB(1, 1.0 - degree, 1.0 - degree)
                elif name.startswith('Forebrain'):
                    cnv.setFillColorRGB(1.0 - degree, 1.0 - degree, 1)
                else:
                    cnv.setFillColorRGB(1.0 - degree, 1.0 - degree, 1.0 - degree)
                cnv.rect(x0, y, x1 - x0, height, 0, 1)
                pass
            y -= height
        #
        x += width + spacer
    cnv.save()
#        xcnv = lambda x: float(x - start) / (stop - start) * width + gx


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', nargs='+', default=None, metavar='name1 name2 ...', help='titles')
    parser.add_argument('-l', nargs='+', metavar='chr*:start-stop chr*:start:stop...', help='genome location')
    parser.add_argument('-o', default='profile.pdf', help='pdf filename', metavar='filename')
    parser.add_argument('-B', default=None, help='background', metavar='name')
    parser.add_argument('-y', default=None, metavar='min,max')
    parser.add_argument('-g', default='200,100', metavar='geometry (width, height[,spacer])')
    parser.add_argument('-G', default=None, metavar='filename', help='GTF')
    parser.add_argument('-w', type=int, default=250, metavar='number', help='window size')
    parser.add_argument('-d', default=None, metavar='directory', help='data base directory')
    parser.add_argument('-s', type=int, default=0, help='smooth level', metavar='number')
    parser.add_argument('-b', help='BED file', default=None, nargs='+', metavar='bed file')
    parser.add_argument('--bar', action='store_true', help='line or bar')
    parser.add_argument('--verbose', action='store_true')
    #parser.add_argument('--barcode', action='store_true')
    parser.add_argument('--barcode', default=None, help='show as barcode when RGB is given')
    parser.add_argument('--overlay', action='store_true', help='overlay graphs')
    parser.add_argument('--log', action='store_true', help='log ratio with -B option')
    parser.add_argument('--with-ticker', action='store_true')
    parser.add_argument('--vspace', type=float, default=10.0)
    parser.add_argument('--hspace', type=float, default=70.0)
    #parser.add_argument('--stat', action='store_true')
    parser.add_argument('--chipchip', action='store_true', help='draw ChIP-chip data in the same format')
    parser.add_argument('-c', nargs='+', default=None, help='ChIP-chip data')
    args = parser.parse_args()

#python ~/Research/Polycomb/python/chipseqgraph.py -d atacmm10 -b ATAC/cTreg_over_eTreg.bed ATAC/eTreg_over_cTreg.bed -n Tn_ cTreg_ Teff_ eTreg_ -o hoge.pdf -G /Data/Genomes/M_musculus/mm10/genes.gtf  -y 0,10 -s 1 -g 150,20 --vspace 0 -l chr17:65920000-65980000 chr1:164040000-164090000 chr2:11620000-11700000
#
#
#
#
#
#
#


    gtf = None
    overlay = args.overlay
    bedfiles = args.b
    barheight = 8
    markedregions = []
    background = args.B
    mode_log = args.log
    use_ticker = args.with_ticker
    vspacer = args.vspace
#    if args.b:
#        markedregions = load_bed(args.b)

    if args.barcode is not None:
        if re.match('\\d+,\\d+,\\d+$', args.barcode):
            color = reportlab.lib.colors.toColor('rgb({})'.format(args.barcode))
        else:
            color = reportlab.lib.colors.toColor(args.barcode)
        barcode_color = [color.red, color.green, color.blue]
    else:
        barcode_color = None
    chipchip = args.chipchip
    if chipchip:
        if args.c is None:
            raise Exception('chipchip require input files')
        else:
            samples = args.c
    else:
        samples = args.n
    names = args.n
    if names is None:
        if chipchip is False:
            raise Exception('no name')
        names = args.c

    # if args.barcode:
    #     draw_barcode_graphs(args)
    #     exit()

    geo = [float(x_) for x_ in args.g.split(',')]
    width = geo[0]
    height = geo[1]
    spacer = args.hspace#70
    division = args.w
    if len(geo) > 2: spacer = geo[2]
    N = (len(args.n) if args.n is not None else 0) + (len(args.c) if args.c is not None else 0)
    num_graphs = len(args.l) * N
    pw = (spacer + width) * (len(args.l) + 1) + spacer
    ph = (spacer + height) * N + spacer * 3
    pagesize = pw, ph
    cnv = reportlab.pdfgen.canvas.Canvas(args.o, pagesize=(pw, ph))
    arrow = cnv.beginPath()
    arrow.moveTo(0,0)
    arrow.lineTo(0, 4)
    arrow.lineTo(10, 4)
    arrow.moveTo(8, 6)
    arrow.lineTo(10, 4)
    arrow.lineTo(8,2)
    colors = ['ReportLabFidRed', 'ReportLabFidBlue', 'ReportLabGreen',
              'orange', 'purple', 'pink', 'olive', 'cadetblue', 'brown', 'navy',
              'gray', 'magenta']
    cnv.setFont('Helvetica', 10)
    if args.y is not None:
        ymin, ymax = [float(x_) for x_ in args.y.split(',')]
    gx = spacer
    total_counts = {}
    if args.n:
        for name in args.n:
            total_counts[name] = count_total_sequences(args.d, name)
            #print(name, total_counts[name])
    for location in args.l:
        #print(location)
        m = re.match('((chr)?[^:]+):((\\+|\\-):)?(\\d+)\\-(\\d+)$', location)
        if m:
            chrom = m.group(1)
            start = int(m.group(5))
            stop = int(m.group(6))
        else:
            raise Exception('cannot resolve position {}'.format(location))
        if start > stop: start, stop = stop, start
#        print(chrom, start, stop)
#        chrom, pos = location.split(':')
#        start, stop = [int(x_.replace(',', '')) for x_ in pos.split('-')]
        xcnv = lambda x: float(x - start) / (stop - start) * width + gx
        if background is not None and not chipchip:
            total_denom = count_total_sequences(args.d, background)
            denominator = convert_dict_to_profile(get_profile(args.d, background, chrom, start, stop), start, stop, division, sum(total_denom.values()))
        else:
            denominator = None
        gy = (height + spacer) * N - height + spacer
        for index, name in enumerate(samples):
            if name not in total_counts or total_counts[name] is None:
                raise Exception('{} is not available'.format(name))
                continue
            if chipchip:
                prof = get_chipchip_profile(name, chrom, start, stop)
                values = convert_dict_to_chipchip_profile(prof, start, stop, division)
            else:
                tc = sum(total_counts[name].values())
                prof = get_profile(args.d, name, chrom, start, stop)
                values = convert_dict_to_profile(prof, start, stop, division, tc)
            #print(index, name, len(values), max(values))
            smin = 0#min(counts.keys())
            smax = (stop - start) // division + 1
            if denominator is not None:
                ratios = []
                for i in range(min(len(values), len(denominator))):
                    v1 = values[i]
                    v2 = denominator[i]
                    if mode_log:
                        if v2 > 0 and v1 > 0:
                            logratio = math.log(v1 / v2) / math.log(2)
                        else:
                            logratio = 0.0
                    else:
                        logratio = v1 - v2
                    ratios.append(logratio)
                    pass
                values = ratios
            if len(values) <= 1:
                sys.stderr.write('too few elements {}\n'.format(len(values)))
                continue
            if args.s > 0:
                values = smoothen(values, args.s)
            if args.y is None:
                ymin = 0
                ymax = max(values) * 1.1
            #print(ymin, ymax)
            ycnv = lambda y: float(max(ymin, min(ymax, y)) - ymin) / (ymax - ymin) * height + gy
            path = None
            ypoints = []
            if barcode_color:
                bw = xcnv(division) - xcnv(0)
                vx = xcnv(start)
                for j, val in enumerate(values):
                    vx = xcnv(j * division + start)
                    degree = (val - ymin) / (ymax - ymin)
                    c_ = [(1 - degree * (1 - barcode_color[i])) for i in range(3)]
                    cnv.setFillColorRGB(c_[0], c_[1], c_[2])
                    cnv.rect(vx, gy, bw, height, 0, 1)
            else:
                for j, val in enumerate(values):
                    vx = xcnv(j * division + start)
                    if smin <= j < smax:
                        vy = ycnv(val)
                        ypoints.append([vx, vy])
                        #print(j, val, vx, vy)
                        if path is None:
                            path = cnv.beginPath()
                            path.moveTo(vx, gy)
                            pass
                        path.lineTo(vx, vy)
                if args.bar:
                    bw = xcnv(division) - xcnv(0)
                    cnv.setFillColor(colors[index % len(colors)])
                    for vx, vy in ypoints:
                        cnv.rect(vx, gy, bw, vy - gy, 0, 1)
                elif path is not None:
                    path.lineTo(vx, gy)
                    path.close()
                    cnv.setStrokeColor(colors[index % len(colors)])
                    cnv.drawPath(path)
            if not overlay or index == 0:
                cnv.setStrokeColor('black')
                cnv.rect(gx, gy, width, height)
                #ticks
                #print(use_ticker)
                if use_ticker:
                    ticks = cnv.beginPath()
                    ticks.moveTo(gx, gy)
                    xstep = 5000
                    while stop - start > xstep * 10:
                        xstep *= 2
                    x0 = (start // xstep) * xstep
                    x1 = (stop // xstep) * xstep
                    x = x0
                    while x <= x1:
                        if x >= x0:
                            ticks.moveTo(xcnv(x), gy)
                            ticks.lineTo(xcnv(x), gy - 3)
                        x += xstep
                    y = ymin
                    ystep = get_step(ymin, ymax)
                    while y <= ymax:
                        ticks.moveTo(gx, ycnv(y))
                        ticks.lineTo(gx - 3, ycnv(y))
                        y += ystep
                    cnv.drawPath(ticks)
                if index == len(samples) - 1:
                    cnv.setFillColorRGB(.02, .01, .007)
                    cnv.drawString(gx - 20, gy, '{}'.format(ymin))
                    cnv.drawString(gx - 20, gy + height, '{}'.format(ymax))
                    cnv.setFillColorRGB(.005, .01, .007)
                    cnv.drawString(gx - 30, gy - 15, chrom)
                    cnv.drawString(gx, gy - 15, '{}'.format(start))
                    cnv.drawString(gx + width, gy - 15, '{}'.format(stop))
            if not overlay:
                gy -= height + vspacer

        #
        if args.G is not None:
            goffset = 30
            genes = load_genes(args.G, chrom, start, stop)
            left_limit = xcnv(start)
            right_limit = xcnv(stop)
            y0 = gy + height - goffset#- vspacer * 2#+ height + vspacer * 4#+ height - vspacer * .25
            h0 = 4
            shared = []
            for i, gene in enumerate(genes):
                emin = min([x_[0] for x_ in gene.exons])
                emax = max([x_[1] for x_ in gene.exons])
                x0 = xcnv(emin)
                x1 = xcnv(emax)
                layer = 0
                for l, x0_, x1_ in shared:
                    if x0 <= x1_ and x1 >= x0_:
                        layer = l + 1
                shared.append([layer, x0, x1])
                #y0 = gy + height + spacer * .5 - h0 * 2 * layer
                y0 = gy + height - goffset - h0 * 2 * layer
                cnv.setStrokeColorRGB(.10, 0.60, .04)
                cnv.setFillColorRGB(.10, 0.60, .04)
                cnv.line(max(left_limit, x0), y0, min(right_limit, x1), y0)
                for es, ee in gene.exons:
                    x0 = xcnv(es)
                    x1 = xcnv(ee)
                    if x1 >= left_limit and x0 <= right_limit:
                        cnv.rect(x0, y0 - h0 * .5, x1 - x0, h0, 0, 1)
                    pass
                cnv.setFillColorRGB(.5, .7, .5)
                cnv.drawString(x1 + 4, y0, gene.name)

                if gene.orientation == '+':
                    cnv.translate(xcnv(gene.start), y0 + h0)
                else:
                    cnv.translate(xcnv(gene.stop), y0 + h0)
                    cnv.scale(-1,1)
                cnv.drawPath(arrow)
                cnv.resetTransforms()
        if bedfiles is not None:
            boffset = 100
            left_limit = xcnv(start)
            right_limit = xcnv(stop)
            bedcolors = 'violet', 'salmon', 'brown', 'cyan', 'gold', 'pink',
            tickwidth = 0
            # triangle = cnv.beginPath()
            # triangle.moveTo(0, 0)
            # triangle.lineTo(barheight * .7, barheight)
            # triangle.lineTo(barheight * -.7, barheight)
            # triangle.lineTo(0, 0)
            # triangle.close()

            for i, bedfile in enumerate(bedfiles):

                #y0 = gy + height - vspacer * .50 - i * (barheight * 1.5)
                y0 = gy + height - boffset - i * (barheight * 1.5)#height - vspacer * .50 - i * (barheight * 1.5)
                #cnv.setFillColorRGB(.1, .8, .3)
                n_ = 0
                cnv.setFillColor(bedcolors[i % len(bedcolors)])

                for section in load_bed_sections(bedfile, chrom, start, stop):
                    #print(section)
                    #if section.chromosome != chrom: continue
                    x0 = xcnv(section.start)
                    x1 = xcnv(section.end)
                    #print(x0, x1, y0)
                    if x0 <= right_limit and x1 >= left_limit:
                        # if x1 - x0 < 10:
                        #     cnv.translate((x0 + x1) * .5, y0)
                        #     cnv.drawPath(triangle, 0, 1)
                        #     cnv.resetTransforms()
                        # else:
                        tri = cnv.beginPath()
                        tri.moveTo(x0, y0)
                        tri.lineTo(x1, y0)
                        tri.lineTo(x1 + tickwidth, y0 + barheight)
                        tri.lineTo(x0 - tickwidth, y0 + barheight)
                        tri.lineTo(x0, y0)
                        tri.close()
                        cnv.drawPath(tri, 0, 1)
#                        cnv.rect(x0, y0, x1 - x0, barheight, 0, 1)
                        n_ += 1
#                print(bedfile, y0, n_)

        gx += spacer + width
        #print(gx, gy)
        pass
    cnv.save()
