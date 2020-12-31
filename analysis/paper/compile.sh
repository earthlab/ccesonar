#!/bin/bash

# compile to pdf
pdflatex main_src
bibtex main_src
pdflatex main_src
pdflatex main_src

# convert pdf images to eps
inkscape input.pdf --export-eps=output.eps