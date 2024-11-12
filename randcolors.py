import random as r
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--start",type=int,default=10)
args = parser.parse_args()
colors = [(r.random(),r.random(),r.random()) for _ in range(90)]

for i,c in enumerate(colors,start=args.start):
	x,y,z=c
	print("\definecolor{{gene%i}}{{rgb}}{{{:.2},{:.2},{:.2}}}".format(x,y,z)%i)
