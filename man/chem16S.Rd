\encoding{UTF-8}
\name{chem16S}
\alias{chem16S}
\alias{chem16S_1}
\title{Community-level chemical metrics}
\description{
Plots from the paper by Dick and Kang (2023).
}

\usage{
  chem16S_1(pdf = FALSE)
}

\arguments{
  \item{pdf}{logical, make a PDF file?}
}

\details{
\code{chem16S_1} makes Figure 1 in the manuscript.
}

\references{
Dick JM and Kang X (2023) \emph{chem16S}: community-level chemical metrics for exploring genomic adaptation to environments. \emph{Bioinformatics} \bold{39}, btad564. \doi{10.1093/bioinformatics/btad564}
}

\examples{
# The function attaches the phyloseq package (with require()),
# so we mark it dontrun to save time in tests
\dontrun{
chem16S_1()
}
}
