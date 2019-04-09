#!/usr/bin/python
#-*-coding:utf-8-*-

"""Sequence logo
"""

import os, sys, re, math

class PWM( object ):
    "Positional weight matrix providing index proposed by Pan and Pan (2007)"
    def __init__( self, motif_id, name, species, accession=None, medline=None ):
        self.__id = motif_id
        self.__name = name
        self.__accession = accession
        self.__medline = medline
        self.__species = species
        self.__matrix = None
        self.__normalized = None
        self.__max_frequencies = None
        self.__length = -1
        self.__vmax   = 0.0
        pass

    def generate_reversed(self):
        obj = PWM(self.__id, self.__name + '(reversed)', self.__species, self.__accession, self.__medline)
        obj.__vmax = self.__vmax
        obj.__length = self.__length
        rev = [None] * 4
        for row in range(4):
            r = self.__matrix[3 - row]
            rev[row] = [v for v in reversed(r)]
            pass
        obj.set_matrix(rev)
        return obj
    
    def set_matrix( self, value ):
        """ 1,1,1,1/2,2,2,3/ """
        if isinstance(value, str):#type( value ) in types.StringTypes:
            matrix = []
            rows = value.split( '/' )
            for row in rows:
                columns = [ int( x ) for x in row.split( ',' ) ]
                matrix.append( columns )
                pass
            value = matrix
            pass
        self.__matrix = value
        self.__normalized = [[],[],[],[]]
        self.__max_frequencies = []
        self.__length = len( value[ 0 ] )
        for pos in range( 0, self.__length ):
            count = 0
            for i in range( 4 ): 
                count += value[ i ][ pos ]
                pass
            max_freq = 0.0
            for i in range( 4 ): 
                freq = float( value[ i ][ pos ] ) / count 
                if freq > max_freq: max_freq = freq
                self.__normalized[ i ].append( freq )
                pass
            self.__max_frequencies.append( max_freq )
            pass
        pass
    matrix = property(lambda s:s.__matrix)

    def evaluate_weight_sequence(self, sequence,shift=0,  background=None, minspan=6):
        """Motified Pan's method for PWM vs PWM
        おおよそ0.4より大きいとモチーフだと判定できるらしい
        """
        if background is None:
            background = [ 0.25 ] * 4
            pass
        invb =[  float( sum( background ) ) / pb for pb in background ]
        vmax = 0.0
        maxratio = 0.0
        for pos in range(len(sequence)):
            score = 0.0
            lr1 = lr2 = 0.0
            span = min(len(sequence) - shift, self.__length - pos)
            if span < minspan:
                break
            bestscore = 0
            for i in range(span):
                n = sum(sequence[i + shift])
                r1 = 0.0
                r2 = 0.0
                for j in range(4):
                    weight = sequence[i + shift][j] * invb[j]
                    pm = self.__normalized[j][i + pos]
                    r1 += pm * weight
                    r2 += weight
                r2 *= self.__max_frequencies[pos + i]
                # print('{}\t{:.2f}\t{:.2f}'.format(pos + i, r1, r2))                
                if r1 <= 0.0 or r2 <= 0.0:
                    score = 0
                    break
                score += math.log(r1)
                bestscore += math.log(r2)
            ratio = score / bestscore
            # print('{}:{:.4f}'.format(pos, ratio))
            if ratio > maxratio:
                maxratio = ratio

        return maxratio

    def evaluate_sequence( self, sequence, start=0, background=None ):
        """Pan & Pan (2007)による評価法
        おおよそ0.4より大きいとモチーフだと判定できるらしい
        """
        if background is None:
            background = [ 0.25 ] * 4
            pass
        invb =[  float( sum( background ) ) / pb for pb in background ]
        pos = 0
        vmax = 0.0
        v = 0.0
        print(self.__max_frequencies, self.__normalized)
        while pos < self.__length:
            base = sequence[ start + pos ]
            if base == 'A':
                n = 0
            elif base == 'C':
                n = 1
            elif base == 'G':
                n = 2
            elif base == 'T':
                n = 3
            else:
                return 0.0
                pass
            pm = self.__normalized[ n ][ pos ]
            if pm <= 0.0: return 0.0
            r1 = math.log( pm * invb[ n ] )
            r2 = math.log( self.__max_frequencies[ pos ] * invb[ n ] )
            print('{}\t{:.2f}\t{:.2f}\t{:.2f}'.format(pos, pm, r1, r2))
            v += r1
            vmax += r2
            pos += 1
            pass
        print(v, vmax)
        if v < 0.0:
            index = 0.0
        else:
            index = v / vmax
            pass
        return index

    name    = property( lambda s: s.__name )
    id      = property( lambda s: s.__id )
    species = property( lambda s: s.__species )
    length  = property( lambda s: s.__length )
    matrix  = property( lambda s: s.__matrix )

    def __repr__( self ):
        output = "%s\t%s\t%s\n" % ( self.id, self.name, self.species )
        for i, row in enumerate( self.__matrix ):
            line = '{} |'.format('ACGT'[i])
            # line = "%s |" % ( "ACGT"[ i ] )
            for datum in row:
                line += '{}'.format(datum)#"%4d" % ( datum )
                pass
            output += line + '\n'
            pass
        return output

    @classmethod
    def load_motifs( self, species = 'Mus musculus' ):
        import pysqlite2.dbapi2 as sqlite
        cnx = sqlite.connect( '/Users/takaho/Research/TFMotif/db/jaspar.db' )
        cur = cnx.cursor()
        if species is None:
            sqlcom = "select id, name, species, accession, medline, matrix from JASPAR"
        else:
            sqlcom = "select id, name, species, accession, medline, matrix from JASPAR where species=\"%s\"" % ( species )
            pass
        cur.execute( sqlcom )
        motifs = {}
        for motif_id, name, species, accession, medline, matrix in cur.fetchall():
            motif = PWM( motif_id, name, species, accession, medline )
            motif.set_matrix( matrix )
            motifs[ motif.id ] = motif
            #motifs.append( motif )
            #print motif
            pass
        cur.close()
        cnx.close()
        return motifs

    @classmethod
    def get_preset( cls, name ):
        label = name.lower()
        if label == 'myc':
            # JASPER database
            motif = pwm.PWM('MA0147.1', 'Myc', 'Mus musculus')
            motif.set_matrix("67,34,8,217,0,19,9,0,0,45/96,53,219,4,212,2,44,2,17,107/36,130,0,5,3,204,0,216,183,24/28,10,0,1,12,2,174,9,27,51")
            return motif
        elif label == 'oct4' or label == 'pou5f1':
            # Jin et al., Genome Research, 2007
            #   doi: 10.1101/gr.6006107 
            
            motif = PWM( "Jin2007", "Oct4", "Homo sapiens", "---", "10.1101/gr.6006107" )
            matrix = ( ( 331, 881, 92, 53, 100, 954, 833, 363, 301, 321 ),
                       ( 247, 59, 0, 50, 730, 47, 44, 173, 231, 201 ),
                       ( 217, 0, 64, 839, 110, 0, 57, 237, 198, 209 ),
                       ( 207, 62, 846, 60, 62, 1, 68, 229, 272, 271 ) )
            motif.set_matrix( matrix )
            return motif
        elif label == 'znf217':
            # Krig et al., 2007 The journal of Biological Chemistry
            motif = PWM( 'znf', 'Znf217', 'Homo sapiens', '---', '10.1074/jbcM611752200' )
            matrix = ( ( 69, 6, 2, 0, 11, 31, 63, 6 ),
                       ( 15, 0, 0, 97, 67, 18, 15, 88 ),
                       ( 1, 0, 0, 0, 0, 32, 13, 3 ),
                       ( 12, 91, 95, 0, 19, 16, 6, 0 ) )
            motif.set_matrix( matrix )
            return motif
        elif label == 'p53' or label == 'tp53':
            motif = PWM( "Funk1992", "TP53", "Homo sapiens", "---", "medline:1588974" )
            matrix = ((5, 3, 4, 5,13, 0,17, 0, 0, 0, 0, 0, 1, 1, 4, 1,15, 2, 1, 1),
                      (8, 7, 0, 0, 0,17, 0, 0, 0,11,16,16, 0, 0, 0,14, 0, 0, 1, 2),
                      (2, 6,13,12, 4, 0, 0, 0,17, 0, 0, 0,15,14,13, 2, 0, 0,14, 1),
                      (2, 1, 0, 0, 0, 0, 0,17, 0, 6, 1, 1, 1, 2, 0, 0, 2,15, 1,13))
            motif.set_matrix( matrix )
            return motif
        raise Exception( "no data for %s" % name  )

    @classmethod
    def generate_database(cls, directory = 'data/JASPAR_CORE_2008' ):
        "JASPAR_COREのデータを読み取ってデータベースに登録する"
        import pysqlite2.dbapi2 as sqlite
        cnx = sqlite.connect( 'db/jaspar.db' )
        cur = cnx.cursor()
        try:
            cur.execute( 'drop table JASPAR' )
        except:
            pass
        cur.execute( "create table if not exists JASPAR ( id not null primary key, name, species, accession, medline, matrix ) " )
        cur.execute("delete from JASPER")
        filename_list = os.path.join( directory, 'matrix_list.txt' )
        with open(filename_list) as fi:
            for line in fi:#open( filename_list ):
                items = line.strip().split( '\t' )
            #print items
                motif_id = items[ 0 ]
                name = items[ 2 ]
                species = None
                accession = None
                medline = None
                for data in items[ 4 ].split( ';' ):
                    try:
                        label, value = data.strip().split( ' ', 1 )
                        if label == 'species':
                            species = value.strip( '"' )
                        elif label == 'acc':
                            accession = value.strip( '"' )
                        elif label == 'medline':
                            medline = value.strip( '"' )
                            pass
                        pass
                    except:
                        pass
                    pass
                if species is None or medline is None:
                    continue
                sqlcom = "insert into JASPAR ( id, name, species, accession, medline ) values( \"%s\", \"%s\", \"%s\", \"%s\", \"%s\" )" % ( motif_id, name, species, accession, medline )
            #print sqlcom
                cur.execute( sqlcom )
                pass
            pass
        cur.close()

        cur.execute( "select id, name from JASPAR" )
        for motif_id, name in cur.fetchall():
            filename = os.path.join( directory, motif_id + '.pfm' )
            #print filename
            with open(filename) as fi:
                lines = [ re.sub( '\\s+', ',', line.strip() ) for line in fi]#open( filename ).readlines() ]
                pass
            value = lines[ 0 ] + '/' + lines[ 1 ] + '/' + lines[ 2 ] + '/' + lines[ 3 ]
            sqlcom = "update JASPAR set matrix=\"%s\" where id=\"%s\"" % ( value, motif_id )
            #print sqlcom
            cur.execute( sqlcom )
            pass

        cnx.commit()
        cnx.close()
        pass
    pass

class SequenceLogo(object):
    @classmethod
    def calculate_information(cls, nucleotides, base=2.0):
        """Calculate information entropy using Shannon's method """
        import math
        nucleotides = [x_ for x_ in nucleotides if x_ > 0]#filter(lambda x: x >= 0, nucleotides)
        num = sum(nucleotides)
        probs = [float(n_) / num for n_ in nucleotides]
        p = 0
        # print('probs : ', probs, nucleotides, num)
        for pr in probs:
            if pr == 0.0: continue
            p -= math.log( pr ) * pr
            pass
        return p / math.log(base)

    @classmethod
    def draw_on_pdf(cls, canvas, frequencies, position=None, size=None, colors=None, fontname='Helvetica'):
        """Draw sequence logo on pdf canvas.
        canvas   : PDF canvas provided by Reportlab PDF library
        position : tuple of (X, Y) to be drawn (default is (100,300))
        size     : tuple of dimension of logo (default is (N * 10, 100), N = length)
        colors   : dictionary of base to color name
                   default  A:brown, C:navy, G:olive, T:green
        fontname : name of font (default is Helvetica)               
        """
        N = len(frequencies)
        if colors is None:
            colors = { 'A':'brown', 'C':'navy', 'G':'olive', 'T':'green'}
        if position is None: position = (100,300)
        if size is None: size = (N * 10, 100)

        from reportlab.pdfbase import pdfmetrics
        from reportlab.pdfbase.pdfdoc import PDFError, PDFDocument
        #from reportlab.pdfbase.ttfonts import TTFont

        font_ac, font_dec = pdfmetrics.getAscentDescent(fontname)

        def sort_nucleotides(nucleotides):
            nucs = { 'A':nucleotides[0],
                     'C':nucleotides[1],
                     'G':nucleotides[2],
                     'T':nucleotides[3] }
            sorted_nucs = []
            for base in sorted(nucs.keys(), key=lambda a:nucs[a]):#lambda a, b: -(nucs[a] > nucs[b])):
                sorted_nucs.append((base, nucs[base]))
                pass
            return sorted_nucs

        canvas.saveState()
        fontsize = 1.0
        canvas.setFont('Helvetica', fontsize)
        dx = float(size[0]) / float(N)
        font_coeff = 1.0 / (font_ac * 0.001)
        scale_x = dx 
        for i, nucleotides in enumerate(frequencies):
            # print(i, cls.calculate_information(nucleotides))
            height_ratio = (2.0 - cls.calculate_information(nucleotides)) / 2.0 * size[1]\
                / fontsize
            x = i * dx + position[0]
            sorted_nucs = sort_nucleotides(nucleotides)
            total_count = sum(nucleotides)
            y = position[1]
            for base, count in sorted_nucs:
                if count == 0: continue
                canvas.saveState()
                canvas.translate(x, y)
                canvas.setFont(fontname, fontsize)
                scale_y = float(count) / total_count * height_ratio
                # print('{}:{}\tyscale  = {:.2f}, {:.2f}'.format(base, count, scale_y, height_ratio))
                canvas.scale(scale_x * 1.5, scale_y * font_coeff)
                canvas.setFillColor(colors[base])
                canvas.drawString(canvas.stringWidth(base) * -0.5, 0, base)
                canvas.restoreState()
                y += scale_y * fontsize
                pass
            pass
        canvas.restoreState()
        #canvas.rect(position[0], position[1], size[0], size[1])
        pass
    pass

if __name__ == '__main__':
    pwm = PWM('test', 'test', 'at')
    pwm.set_matrix('1,1,0,0,1/0,1,1,0,0/0,1,9,0,2/0,0,0,1,1')
    print(pwm)
    print(pwm.evaluate_sequence('ACGTT'))
    mat = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[0,1,0,1]]#[1,1,1,1], [0,1,0,0], [0,0,1,0], [1,0,0,1], [1,1,0,0], [0,1,0,0]]
    score = pwm.evaluate_weight_sequence(mat)#sequence)
    print(score)

    exit()        
    import random
    from reportlab.pdfgen.canvas import Canvas
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.pdfdoc import PDFError, PDFDocument
    from reportlab.pdfbase.ttfonts import TTFont

#     c = Canvas('logo.pdf')
#     fontname = 'Helvetica'
#     font = pdfmetrics.getFont(fontname)#'Helvetica')#Font('Helvetica')
#     print font
#     print pdfmetrics.getAscent(fontname)
#     print pdfmetrics.getDescent(fontname)
#     print pdfmetrics.getAscentDescent(fontname)

#     font = pdfmetrics.getFont('Courier')
    
#     print pdfmetrics.getRegisteredFontNames()
#     fontsize = 25
#     c.setFont(fontname, fontsize)
#     text = "Font test FONT TEST"
#     c.drawString(100, 300, text)
#     c.rect(100, 300, font.stringWidth(text, fontsize), pdfmetrics.getAscent(fontname) * 0.001 * fontsize )
#     c.save()

#     exit()

    import re, sys, os
    info = None
    alias = {'A':0, 'C':1, 'G':2, 'T':3}
    info = None
    counts = None
    c = Canvas('logo.pdf')
    x = 50
    y = 700
    with open(sys.argv[1]) as fi:
        for line in fi:
            if line.startswith('Score'):
                info = line.strip()
                counts = None
            elif len(line) <= 1:
                info = None
            elif info is not None:
                if line.find('**') >= 0:
                    c.drawString(x, y + 50, info)
                    draw_on_pdf(c, counts, position=(x, y), size=(len(counts) * 30, 50))
                    y -= 75
                    info = None
                    pass
                else:
                    items = re.split('\\s+', line.strip())#.split('\t')
                    print(items)
                    if len(items) > 4:
                        pattern = re.sub('[acgtn]', '', items[2])
                        if counts is None:
                            counts = [ [0] * 4 for i in range(len(pattern))]
                            pass
                        for i in range(len(pattern)):
                            base = pattern[i]
                            if base not in alias: break
                            counts[i][alias[base]] += 1
                            pass
                        pass
                    pass
                pass
            pass
        pass
    c.save()
    exit(0)
    
    c = Canvas('logo.pdf')
    data = ( ( 10, 0, 0, 0 ),
             ( 8, 1, 1, 0 ),
             ( 1, 2, 5, 3 ),
             ( 0, 0, 5, 5 ),
             ( 2, 3, 2, 3 ),
             ( 1, 100, 2, 0 ))
    draw_on_pdf(c, data, fontname='Courier' )
    data = ( (1,2,1,7), (10,0,1,0), (0,0,0,11), (0,0,11,0), (5,4,0,1), )
    draw_on_pdf(c, data, position=(50,200), size=(400,100))
    c.save()

