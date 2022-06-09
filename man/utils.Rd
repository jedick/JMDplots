\encoding{UTF-8}
\name{utils}
\alias{utils}
\alias{addexif}
\alias{addalpha}
\alias{cldfun}
\title{Some Utilities for the Package}
\description{
  Utility functions for the package.
}

\usage{
  addexif(name, title, paperref)
  addalpha(col, alpha)
  cldfun(ZClist, bp, dy)
}

\arguments{
  \item{name}{character, function name (added to Creator in PDF metadata)}
  \item{title}{character, title for PDF}
  \item{paperref}{character, reference (added to Subject in PDF metadata)}
  \item{col}{color specification (see \code{\link{col2rgb}})}
  \item{alpha}{character, hexadecimal value of alpha}
  \item{ZClist}{list, values of \ZC}
  \item{bp}{output of \code{\link{boxplot}}}
  \item{dy}{numeric, how far above box to plot letters}
}

\details{
\code{addexif} is a utility function to write add PDF metadata using the \samp{exiftool} system command.
The file name is constructed from the function name (i.e. \code{name.pdf}).
If \samp{exiftool} is not available, this function does nothing.

\code{addalpha} adds transparency to a color by appending the value of alpha to the hexadecimal string equivalent of the color given in \code{col}.

\code{cldfun} adds significant difference letters to a boxplot.
\code{ZClist} should have the values of \ZC with names for each group.
}
