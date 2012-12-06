# Makefile
#

LATEX = latex
BIBTEX = bibtex
DVIPS = dvips
PS2PDF = ps2pdf14

PDF = $(addsuffix .pdf, $(basename $(shell grep -l '\\begin{document' *.tex)))
BIB = $(addsuffix .bbl, $(basename $(shell grep -l '\\nocite{\|\\cite{' *.tex)))


all: $(BIB) $(PDF)

clean:
	rm -f *.blg *.log *.nav *.out *.snm *.toc *.dvi *.aux

release: all clean

distclean: clean
	rm -f $(BIB) $(PDF)

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
