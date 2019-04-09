#!/usr/bin/python

"""Sequence logo
"""

def calculate_information(nucleotides, base=2.0):
    """Calculate information entropy using Shannon's method """
    import math
    nucleotides = filter(lambda x: x >= 0, nucleotides)
    num = sum(nucleotides)
    probs = [float(n_) / num for n_ in nucleotides]
    p = 0
    for pr in probs:
        if pr == 0.0: continue
        p -= math.log( pr ) * pr
        pass
    return p / math.log(base)
#    return 2.0 - p / math.log( 2.0 )

# def draw_png(frequencies, size=None):
#     import Image, ImageDraw
#     if size is None:
#         height = 200
#         width = len(frequencies) * 40
#     else:
#         width = size[0]
#         height = size[1]
#         pass
#     class ImageObject(object):
#         def __init__(self,img):
#             self.__img = img
#             self.__filename = tempfile.mktemp('.png')
#             fo = open(self.__filename, 'w')
#             self.__img.save(fo)
#             fo.close()
#             pass
#         def __del__(self):
#             if self.__filename is not None:
#                 os.path.unlink(self.__filename)
#                 pass
#             pass
#         filename = property(lambda s: s.__filename)
#         image = property(lambda s: s.__img)
#         pass

#     img = Image.new('RGB', (width, height))
#     for i, nucs in enumerate(frequencies):
#         x0 = i * width / len(frequencies)
#         x1 = 
#         pass
#     seqlogo = ImageObject(img)
#     return seqlogo

def draw_on_pdf(canvas, frequencies, position=None, size=None, colors=None, fontname='Helvetica'):
    """Draw sequence logo on pdf canvas.
    canvas   : PDF canvas provided by Reportlab PDF library
    frequencies : [[a1,c1,g1,t1],[a2,c2,g2,t2],....]
    position : tuple of (X, Y) to be drawn (default is (100,300))
    size     : tuple of dimension of logo (default is (N * 10, 100), N = length)
    colors   : dictionary of base to color name
               default  A:red, C:blue, G:yellow, T:green
    fontname : name of font (default is Helvetica)               
    """
    N = len(frequencies)
    if colors is None:
        colors = { 'A':'crimson', 'C':'ReportLabFidBlue', 
                   'G':'gold',
                   'T':'ReportLabGreen'}
    if position is None: position = (100,300)
    if size is None: size = (N * 10, 75)

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
        for base in sorted(nucs.keys(), lambda a, b: -(nucs[a] > nucs[b])):
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
        height_ratio = (2.0 - calculate_information(nucleotides)) / 2.0 * size[1]\
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

if __name__ == '__main__':
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

