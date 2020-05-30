#!/bin/bash

# compile to pdf
pdflatex main_src
bibtex main_src
pdflatex main_src
pdflatex main_src
