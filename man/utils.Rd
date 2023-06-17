\encoding{UTF-8}
\name{utils}
\alias{utils}
\alias{addexif}
\alias{addalpha}
\alias{cldfun}
\alias{hyphen.in.pdf}
\title{Some Utilities for the Package}
\description{
  Utility functions for the package.
}

\usage{
  addexif(name, title, paperref)
  addalpha(col, alpha)
  cldfun(Zclist, bp, dy)
  hyphen.in.pdf(x)
}

\arguments{
  \item{name}{character, function name (added to Creator in PDF metadata)}
  \item{title}{character, title for PDF}
  \item{paperref}{character, reference (added to Subject in PDF metadata)}
  \item{col}{color specification (see \code{\link{col2rgb}})}
  \item{alpha}{character, hexadecimal value of alpha}
  \item{Zclist}{list, values of \Zc}
  \item{bp}{output of \code{\link{boxplot}}}
  \item{dy}{numeric, how far above box to plot letters}
  \item{x}{character, any text object}
}

\details{
\code{addexif} is a utility function to write add PDF metadata using the \samp{exiftool} system command.
The file name is constructed from the function name (i.e. \code{name.pdf}).
If \samp{exiftool} is not available, this function does nothing.

\code{addalpha} adds transparency to a color by appending the value of \code{alpha} to the hexadecimal representation of the color given in \code{col}.

\code{cldfun} adds significant difference letters to a boxplot.
\code{Zclist} should have the values of \Zc with names for each group.

If a PDF device is active, \code{hyphen.in.pdf} replaces all instances of ASCII hyphen (\verb{-}) in \code{x} with the Unicode hyphen (\verb{\uad}).
This is done so that hyphens don't get displayed as minus signs in R graphics output, but it only works for the PDF device.
}
