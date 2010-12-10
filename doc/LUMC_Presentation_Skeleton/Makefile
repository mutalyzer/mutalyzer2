# General Makefile for LaTeX documents by J. F. J. Laros.
# Last alteration on 15-08-2010.
#
# The packages texlive-base-bin, texlive-latex-base and ghostscript should
#   be installed.
#

LATEX = latex
BIBTEX = bibtex
DVIPS = dvips
PS2PDF = ps2pdf14

SRC := $(shell grep -H '\\begin{document' *.tex | cut -f 1 -d '.')

all: $(SRC)

$(SRC): $(SRC).tex
	$(LATEX) $^
	$(LATEX) $^ 
	$(DVIPS) $(SRC).dvi -o $(SRC).ps

bib: $(SRC).tex
	$(LATEX) $^
	$(BIBTEX) $(SRC)

release: $(SRC) clean
	$(PS2PDF) $(SRC).ps
	rm -f $(SRC).ps

clean:
	rm -f *.aux $(SRC).blg $(SRC).dvi $(SRC).log

distclean: clean
	rm -f $(SRC).ps $(SRC).pdf
