#!/bin/bash

latex paper.tex
bibtex paper.aux
latex paper.tex
dvipdf paper.dvi
evince paper.pdf
