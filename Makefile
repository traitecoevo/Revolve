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

test:
	make -C inst/tests

.PHONY: all install doc clean check
