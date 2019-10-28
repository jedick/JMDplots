\encoding{UTF-8}
\name{utils}
\alias{utils}
\alias{addexif}
\title{Some Utilities for the Package}
\description{
  Utility functions for the package.
}

\usage{
  addexif(name, title, paperref)
}

\arguments{
  \item{name}{character, function name (added to Creator in PDF metadata)}
  \item{title}{character, title for PDF}
  \item{paperref}{character, reference (added to Subject in PDF metadata)}
}

\details{
\code{addexif} is a utility function to write add PDF metadata using the \samp{exiftool} system command.
The file name is constructed from the function name (i.e. \code{name.pdf}).
If \samp{exiftool} is not available, this function does nothing.
}
