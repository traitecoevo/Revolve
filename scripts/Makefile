## Rscript really is a festering horrible mess for requiring things
## like this.
RSCRIPT_PKGS := $(shell Rscript -e 'writeLines(Sys.getenv("R_DEFAULT_PACKAGES"))')
RSCRIPT = Rscript --default-packages="${RSCRIPT_PKGS},methods"

runall:
	for f in `ls *.R`; do echo $$f; ${RSCRIPT} $$f; done
