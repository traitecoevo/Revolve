\name{rstar_plot}
\alias{rstar_plot}
\alias{rstar_trajectory}
\title{Plot Rstar ZNGIs}
\usage{
rstar_plot(m, x, xlim = c(0, 1), ylim = c(0, 1), col = seq_len(ncol(x)))

rstar_trajectory(m, x, y = rep(1, ncol(x)), col = seq_len(ncol(x)),
  eps = 1e-06, t_max = 100, t_len = 101, S = NULL, col_died = "grey")
}
\arguments{
  \item{m}{Rstar model, made by \code{make_rstar}}

  \item{x}{Matrix of parameters indicating species traits}

  \item{xlim}{X axis limits (resource 1)}

  \item{ylim}{Y axis limits (resource 2)}

  \item{col}{Vector of colours along species}

  \item{y}{Densities of the species}

  \item{eps}{Density at which species are considered
  extinct}

  \item{t_max}{Maximum time to run simulation until}

  \item{t_len}{Number of time intervals to run simulation
  over}

  \item{S}{Resource supply rates (and initial resource
  levels)}

  \item{col_died}{Colour of communities that go extinct.}
}
\description{
Plot Rstar ZNGIs
}
\details{
Will change
}
\author{
Rich FitzJohn
}

