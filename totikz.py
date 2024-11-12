from argparse import ArgumentParser
from math import pi

CHR_CIRCULAR = ')'
CHR_LINEAR = '|'
ORIENT_POSITIVE = '+'
ORIENT_NEGATIVE = '-'

COLOR_MAP = dict(zip([str(i) for i in range(1,1000)], ['gene%d'%i for i in range(1,100)] ))

SUBID_SEP = '.'
dd = lambda x: '' if x==ORIENT_POSITIVE else x
EXTREMITY_TAIL = 't'
EXTREMITY_HEAD = 'h'

class Marker:
    def __init__(self, gene, mid, direction=ORIENT_POSITIVE):
        self.gene = gene
        self.mid = mid
        self.direction = direction
    def __str__(self):
        return dd(self.direction) +str(self.gene)
    def get_extremities(self,directionless=False):
        exts = [((self.mid),EXTREMITY_TAIL), ((self.mid),EXTREMITY_HEAD)]
        if self.direction == ORIENT_POSITIVE or directionless:
            return exts
        else:
            return exts[::-1]

def makemarker(s, fams):
    '''
    Create a marker from an unimog-marker entry, while extending the family dictionary.
    '''
    s = s.strip()
    orient = ORIENT_POSITIVE
    if s.startswith(ORIENT_NEGATIVE):
        orient = ORIENT_NEGATIVE
        s = s[1:].strip()
    elif s.startswith(ORIENT_POSITIVE):
        s = s[1:].strip()
    if not s in fams:
        fams[s] = []
    fams[s].append(s)
    return Marker(s,'{}{}{}'.format(s,SUBID_SEP,len(fams[s])),direction = orient)

def readGenomes(fl):
    '''
    Read multiple genomes from an unimog-file.
    '''
    gen_name = ''
    chromosomes = []
    genomes = []
    at_least_one = False
    fams = {}
    with open(fl) as f:
        for line in f:
            line = line.strip()
            if line == '':
                continue
            if line.startswith('>'):
                genomes.append((gen_name, chromosomes))
                gen_name = line[1:].strip()
                chromosomes = []
                at_least_one = True
                continue
            if line[-1] not in [CHR_LINEAR, CHR_CIRCULAR]:
                raise AssertionError('%s not a recognized chromosome type (%s for linear or %s for circular).'%(line[-1], CHR_LINEAR, CHR_CIRCULAR))
            chromosomes.append((line[-1],[makemarker(s,fams) for s in line[:-1].split()]))
    if at_least_one:
        genomes = genomes[1:]
        genomes.append((gen_name, chromosomes))
    return genomes


def space_requirement(chrm,ml=2,xpad=1,ypad=1):
    chrt, markers = chrm
    if chrt == CHR_LINEAR:
        return (len(markers)*ml+xpad,ypad)
    else:
        d = len(markers)*ml/pi
        return (d+xpad,d+ypad)

def display_circular(l,x=0,y=0,names=False, markerlen=2):
	circumference = len(l)*markerlen
	radius = circumference/(2*pi)
	x_ = x
	print('%circular chromosome')
	for i,g in enumerate(l):
		color = COLOR_MAP[g.gene]
		startang = -i/len(l)*360
		endang = -(i+1)/len(l)*360
		if g.direction == ORIENT_NEGATIVE:
			startang,endang = endang,startang
		print('\\draw[gene,{color}] ({x},{y})++({startang:.2f}:{radius:.2f}) arc ({startang:.2f}:{endang:.2f}:{radius:.2f});'.format(color=color,x=x_,y=y,startang=startang,endang=endang,radius=radius))
	
def display_linear(l,x=0,y=0,markerlen=2):
     print("%linear chromosome")
     for i,g in enumerate(l):
        color = COLOR_MAP[g.gene]
        if g.direction==ORIENT_POSITIVE:
            rel_start = (i-len(l)/2)*markerlen
            print("\\draw[gene,{color}] ({x},{y})++({relstart:.2f},0) -- ++({markerlen:.2f},0);".format(color=color,x=x,y=y,relstart=rel_start,markerlen=markerlen))
        else:
            rel_start = (i+1-len(l)/2)*markerlen
            print("\\draw[gene,{color}] ({x},{y})++({relstart:.2f},0) -- ++(-{markerlen:.2f},0);".format(color=color,x=x,y=y,relstart=rel_start,markerlen=markerlen))

def t(n):
    A="abcdefghijklmnopqrstuvwxyz"
    if not n<len(A):
        raise AssertionError("By default only up to {} different chromosomes possible.".format(len(A)))
    return A[n]

def display_chromosomes(chromosomes):
    linears = [c for c in chromosomes if c[0]==CHR_LINEAR]
    circulars = [c for c in chromosomes if c[0]==CHR_CIRCULAR]
    max_width_l = max([space_requirement(c)[0] for c in linears],default=0)
    max_width_c = max([space_requirement(c)[0] for c in circulars],default=0)
    middle_l = max_width_l/2
    middle_c = max_width_l + 0.5*(max_width_c)
    linearys=[]
    coord = 0
    for l in linears:
        s = space_requirement(l)[1]
        coord+=0.5*s
        linearys.append(coord)
        coord+=0.5*s
    circularys=[]
    coord = 0
    for c in circulars:
        s = space_requirement(c)[1]
        coord+=0.5*s
        circularys.append(coord)
        coord+=0.5*s
    print('\\begin{tikzpicture}')
    print('\\tikzmath{')
    print('\\genescale = 0.8;')
    print("\\lm = {:.2f};".format(middle_l))
    print("\\cm = {:.2f};".format(middle_c))
    for i,v in enumerate(linearys):
        print("\\ly{} = {:.2f};".format(t(i),-v))
        print("\\lx{} = {};".format(t(i),"\\lm"))
    for i,v in enumerate(circularys):
        print("\\cy{} = {:.2f};".format(t(i),-v))
        print("\\cx{} = {};".format(t(i),"\\cm"))
    print("}")
    for i, c in enumerate(circulars):
         display_circular(c[1],'\\cx%s'%t(i),'\\cy%s'%t(i))
    for i, l in enumerate(linears):
         display_linear(l[1],'\\lx%s'%t(i),'\\ly%s'%t(i))
    print('\\end{tikzpicture}')


    




def main():
    parser = ArgumentParser()
    parser.add_argument('infile')
    args=parser.parse_args()
    gs = readGenomes(args.infile)
    for g in gs:
        print("%Genome '{}'".format(g[0]))
        display_chromosomes(g[1])
    

if __name__ == '__main__':
    main()
