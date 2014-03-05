all: models.pdf

%.pdf: %.md include.tex
	pandoc $< --csl=the-american-naturalist.csl --toc --filter pandoc-citeproc --bibliography=refs.bib --template=include.tex --variable fontsize=12pt --latex-engine=xelatex -o $@

clean:
	rm -f *.pdf
