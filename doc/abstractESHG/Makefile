# Makefile
#

LATEX = latex
BIBTEX = bibtex
DVIPS = dvips
PS2PDF = ps2pdf14
GNUPLOT = gnuplot

PDF = $(addsuffix .pdf, $(basename $(shell grep -l '\\begin{document' *.tex)))
BIB = $(addsuffix .bbl, $(basename $(shell grep -l '\\nocite{\|\\cite{' *.tex)))
EPS = $(addsuffix .eps, $(basename $(shell ls *.gnp)))


all: $(EPS) $(BIB) $(PDF)

clean:
	rm -f *.blg *.log *.nav *.out *.snm *.toc *.dvi *.aux *.vrb

release: all clean

distclean: clean
	rm -f $(PDF) $(EPS)

%.aux: %.tex
	$(LATEX) $^
	rm $(addsuffix .dvi, $(basename $^))

%.bbl: %.aux
	$(BIBTEX) $(basename $^)

%.dvi: %.tex
	$(LATEX) $^
	$(LATEX) $^
	$(LATEX) $^

%.ps: %.dvi
	$(DVIPS) $^ -o $@

%.pdf: %.ps
	$(PS2PDF) $^

%.eps: %.gnp
	$(GNUPLOT) < $<
