# methyPre
Methylation array preprocessing.

This package is downloadable using the following in R:   
`require(devtools);   install_github("metamaden/methyPre")`

##Purpose
This package contains a single wrapper function that serves as an all-in-one workflow assembly. It takes an experiment summary object as used in minfi (eg. RGChannelSet, MethylSet, GenomicRatioSet, etc.) and applies minfi normalizations, filters, and batch correction. The workflow decision tree, including all arguments and operations, is as follows:

![alt text](https://github.com/metamaden/methyPre/blob/master/methypre_workflow1.jpg "methyPre workflow")

##Citations
For more info, see [the minfi Bioconductor page](http://bioconductor.org/packages/release/bioc/html/minfi.html).  For information on cross-reactive and poor quality probes, see [Chen et al 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3592906/) (HM450 platform), and [Pidsley et al 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1066-1) (EPIC platform).
