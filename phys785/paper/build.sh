#!/bin/bash

latex paper.tex
latex paper.tex
dvipdf paper.dvi
evince paper.pdf
