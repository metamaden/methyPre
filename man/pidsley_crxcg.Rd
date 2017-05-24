\name{pidsley_crxcg}
\alias{pidsley_crxcg}
\docType{data}
\title{
Vector of IDs for CpG probes on the EPIC array platform, which are potentially effected by cross-hybridization and genetic variation.
}
\description{
Pidsley et al 2016
}
\usage{data(pidsley_crxcg)}
\format{
  A vector of CpG probe IDs.
}
\details{
}
\source{
}
\references{
Pidsley, Ruth, Elena Zotenko, Timothy J. Peters, et. al. "Critical evaluation of the Illumina MethylationEPIC BeadChip microarray for whole-genome DNA methylation profiling". Genome Biology 2016 17:208.
}
\examples{
# require(wateRmelon); data(melon)
data(pidsley_crxcg)

dim(melon)

#Features  Samples
#    3363       12

dim(melon[!rownames(melon) \%in\% pidsley.crxcg,])

#Features  Samples
#    3115       12
}
\keyword{datasets}
