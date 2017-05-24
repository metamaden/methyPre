\name{illumina_crxcg}
\alias{illumina_crxcg}
\docType{data}
\title{
Vector of IDs for poor-quality CpG probes on EPIC array, as specified by Illumina.
}
\description{
IDs to poor-quality CpG probes on the Illumina EPIC array platform, as ascertained by Illumina in-house in Spring 2017. It is recommended these be excluded in all EPIC array investigations.
}
\usage{data(illumina_crxcg)}
\format{
  A vector of CpG probe IDs.
}
\details{
}
\source{
}
\references{
}
\examples{
# require(wateRmelon); data(melon)
data(illumina_crxcg)

dim(melon)

#Features  Samples
#    3363       12

dim(melon[!rownames(melon) \%in\% illumina.crxcg,])

#Features  Samples
#    3354       12
}
\keyword{datasets}
