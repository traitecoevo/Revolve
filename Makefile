all:

install:
	R CMD INSTALL .

doc:
	make -C doc

clean:
	make -C src clean

check:
	R CMD build .
	R CMD check `ls -1tr Revolve*gz | tail -n1`

.PHONY: all install doc clean check
