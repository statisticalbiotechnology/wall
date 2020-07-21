main=paper

all: ${main}.pdf supplement.pdf coverletter.pdf in_brief.pdf diff.pdf

%.pdf: %.tex
	pdflatex -shell-escape $*.tex
	if ( grep -q citation $*.aux ) ; then \
		bibtex $* ; \
		pdflatex -shell-escape $*.tex ; \
	fi
	pdflatex -shell-escape $*.tex
	#latex $*.tex
	#dvipdfm $*.dvi

diff.tex: ${main}.tex
#	git show quandenser_ncomm_1:./$? > /tmp/$?
	git show d696d9826871:./$? > /tmp/$?
	latexdiff /tmp/$? $? > diff.tex

%.pdf: %.svg
	inkscape $*.svg --export-filename=$*.pdf

clean:
	rm -f *.aux *.bbl *.blg *.log *~
