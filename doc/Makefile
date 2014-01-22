all: models.pdf

%.pdf: %.md include.tex
	pandoc $< --csl=the-american-naturalist.csl --filter pandoc-citeproc --bibliography=refs.bib --template=include.tex --variable mainfont="Times New Roman" --variable sansfont=Arial --variable fontsize=12pt --latex-engine=xelatex -o $@

clean:
	rm -f *.pdf
