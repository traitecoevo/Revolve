##' Utility for taking care of getting/setting/validating parameters.
##'
##' Eventually I'll add a hook here for checking that parameters are
##' allowable, for setting ranges, etc.
##' @title Utility for Parameter Handling
##' @param defaults List of default parameter values.  Allowable
##' parameter names are also taken from this list
##' @param where Environment into which parameters will be
##' set/retrieved.  If ommited a new environment is generated (so
##' access is only via get/set).
##' @return An object of class \code{parameters}
##' @author Rich FitzJohn
##' @export
make_parameters <- function(defaults, where=new.env()) {
  check.names <- function(x) {
    names <- names(x)
    if (length(x) > 0 &&
        (is.null(names) || any(is.na(names)) || any(names == "") ||
         any(duplicated(names))))
      stop("Names must be present, non-empty and non duplicated")
  }
  check <- function(parameters) {
    check.names(parameters)
    if (!all(names(parameters) %in% names))
      stop("Unknown parameters: ", setdiff(names(parameters), names))
    invisible(TRUE)
  }
  set <- function(parameters=list()) {
    check(parameters)
    for (v in names(parameters))
      assign(v, parameters[[v]], where)
  }
  get <- function()
    structure(lapply(names, base::get, where), names=names)

  check.names(defaults)
  names <- names(defaults)
  set(defaults)
  rm(defaults)

  ret <- list(set=set, get=get, names=function() names)
  class(ret) <- "parameters"
  ret
}

##' @S3method names parameters
names.parameters <- function(x, ...)
  x$names()
##' @S3method names<- parameters
`names<-.parameters` <- function(x, value)
  stop("Cannot set names of a parameters object")
