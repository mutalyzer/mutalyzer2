# Generate PDF from LaTeX and BibTeX source.
#
# There are some hacks in here to work with LaTeX compilation, I'll try
# to document them with comments :)
#
# Martijn Vermaat, martijn@vermaat.name

# Configuration
DOCUMENT = mutalyzer-mobile-2013
PDFLATEX = /usr/bin/pdflatex
PDFLATEXFLAGS = -halt-on-error -interaction errorstopmode
BIBTEX = /usr/bin/bibtex

# Just create the PDF!
all: pdf

# Compile BibTeX source file (run this when citations or .bib file change)
bibtex: $(DOCUMENT).bbl

# Create the files
pdf: $(DOCUMENT).pdf

# For the PDF, we need a .tex LaTeX source. Actually, we also need a .bbl
# BibTeX database, but we don't want to regenerate that on every .tex change.
# We compile the file twice to make sure all references are okay.
%.pdf: %.tex %.bbl
	$(PDFLATEX) $(PDFLATEXFLAGS) $<
	while egrep -q -s 'Rerun (LaTeX|to get cross-references right)' $*.log ;\
	    do \
	        $(PDFLATEX) $(PDFLATEXFLAGS) $< ;\
	    done
	$(PDFLATEX) $(PDFLATEXFLAGS) $< # With ntheorem, another run is needed

# For the .bbl BibTeX database we need the .bib BibTeX source and .tex LaTeX
# source.
# Afterwards, we remove the resulting PDF, because we require some additional
# compilation passes that are done by the %.pdf rule.
%.bbl: %.bib %.tex
	$(PDFLATEX) $(PDFLATEXFLAGS) $*
	$(BIBTEX) $*
	$(PDFLATEX) $(PDFLATEXFLAGS) $*
	while egrep -q -s 'Rerun (LaTeX|to get cross-references right)' $*.log ;\
	    do \
	        $(PDFLATEX) $(PDFLATEXFLAGS) $* ;\
	    done
	/bin/rm $*.pdf

# Please make sure we don't kill any sources...
clean:
	rm -f $(DOCUMENT).aux $(DOCUMENT).pdf $(DOCUMENT).log $(DOCUMENT).toc $(DOCUMENT).out $(DOCUMENT).bbl $(DOCUMENT).blg $(DOCUMENT).nav $(DOCUMENT).snm $(DOCUMENT).vrb
