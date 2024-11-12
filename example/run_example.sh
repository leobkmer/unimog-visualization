#!/bin/bash
python3 ../totikz.py genome1.txt > genome1_tikz.tex
pdflatex example.tex
