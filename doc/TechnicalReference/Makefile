# General Makefile for LaTeX documents by J. F. J. Laros.
# Last alteration on 15-10-2009.
#
# The packages texlive-base-bin, texlive-latex-base and ghostscript should
#   be installed.
#

LATEX = latex
BIBTEX = bibtex
DVIPS = dvips
PS2PDF = ps2pdf14

SRC := $(shell grep -H '\\begin{document' *.tex | cut -f 1 -d '.')

BIB := $(shell grep '\\bibliography{' $(SRC).tex > /dev/null && \
               grep '\\cite{' $(SRC).tex)

all: $(SRC)

$(SRC): $(SRC).tex
	$(LATEX) $^ 
ifdef BIB
	$(BIBTEX) $(SRC)
endif
	$(LATEX) $^
	$(LATEX) $^ 
	$(DVIPS) $(SRC).dvi -o $(SRC).ps

release: $(SRC) clean
	$(PS2PDF) $(SRC).ps
	rm -f $(SRC).ps

clean:
	rm -f *.aux $(SRC).bbl $(SRC).blg $(SRC).dvi $(SRC).log $(SRC).toc

distclean: clean
	rm -f $(SRC).ps $(SRC).pdf
