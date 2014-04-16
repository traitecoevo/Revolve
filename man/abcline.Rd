\name{abcline}
\alias{abcline}
\title{Abline Through Point}
\usage{
abcline(x, y, m, ...)
}
\arguments{
  \item{x}{X coordinate of point to pass through}

  \item{y}{Y coordinate of point to pass through}

  \item{m}{Slope of line}

  \item{...}{Arguments passed to \code{abline} (e.g., col)}
}
\description{
Variant of \code{abline} that puts a line through a point
(x,y)
}
\details{
Like \code{abline}, this does not support recycling (uses
\code{abline} internally).
}
\author{
Rich FitzJohn
}

