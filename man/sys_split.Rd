\name{sys_split}
\alias{sys_split}
\title{Split Community}
\usage{
sys_split(sys, dx)
}
\arguments{
  \item{sys}{Community with components 'x' representing
  species traits.  Must be a single value at the moment
  (i.e., single-resident, single-traits)}

  \item{dx}{Amount to perturb the community.  Either a
  scalar (traits will be symmetrically moved this amount)
  or a vector of length 2.}
}
\description{
Split a community into two parts, naively.
}
\details{
This is just a little utility
}
\author{
Rich FitzJohn
}

