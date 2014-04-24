all:

install:
	R CMD INSTALL .

doc:
	make -C doc

clean:
	make -C doc clean

check:
	R CMD build .
	R CMD check --no-manual `ls -1tr Revolve*gz | tail -n1`
	@rm -f `ls -1tr Revolve*gz | tail -n1`
	@rm -rf Revolve.Rcheck

runall:
	make -C scripts runall

test:
	make -C inst/tests

DEVTOOLS_DOCUMENT=devtools::document(roclets=c('namespace', 'rd'))
document:
	@mkdir -p man
	Rscript -e "library(methods); ${DEVTOOLS_DOCUMENT}"

.PHONY: all install doc clean check
