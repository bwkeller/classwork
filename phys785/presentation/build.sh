#!/bin/bash

latex presentation.tex
latex presentation.tex
dvipdf presentation.dvi
evince presentation.pdf
