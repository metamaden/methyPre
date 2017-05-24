\name{chen_crxcg}
\alias{chen_crxcg}
\docType{data}
\title{
Vector of IDs for cross-reactive or polymorphic CpGs on the Illumina HM450 array.
}
\description{
Listed CpG probe IDs were among the 6-percent of array probes shown in Chen et al 2013 to either be polymorphic or co-hybridize with alternate sequences on the array.
}
\usage{data(chen_crxcg)}
\format{
  A vector of CpG probe IDs.
}
\details{
}
\source{
}
\references{
Chen, Y, Mathieu Lemire, Sanaa Choufani et al. "Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray". Epigenetics 2013 Feb 1; 8(2): 203-209.
}
\examples{
# require(wateRmelon); data(melon)
data(chen_crxcg)

dim(melon)

#Features  Samples
#    3363       12

dim(melon[!rownames(melon) \%in\% chen.crxcg,])

#Features  Samples
#    3250       12

}
\keyword{datasets}
