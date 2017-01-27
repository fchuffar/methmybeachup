# methmybeachup
The __methmybeachup__ package provides a set of functions allowing to easily convolves and vizualise methylome signal across genes and CpG islands


## Data

Data set is extracted from GSE51954 Vandiver AR, Irizarry RA, Hansen KD, Garza LA et al. Age and sun exposure-related widespread genomic blocks of hypomethylation in nonmalignant skin. Genome Biol 2015 Apr 16;16:80. PMID: 25886480


## Installation

To get the current development version from github, you need first to install __devtools__. 
Then, install ``methmybeachup``:

```
install.packages("devtools")
devtools::install_github("fchuffar/methmybeachup")
```


## Example

```{r}
layout(1)
library(methmybeachup)

cols = as.numeric(sunexp_design$sex) * 2
names(cols) = rownames(sunexp_design)
gene = genes[6,]
res1 = analyse_meth(gene, sunexp_data, sunexp_platform, cols, PLOT=TRUE)
legend("topright", col=as.numeric(unique(sunexp_design$sex)) * 2, legend=unique(sunexp_design$sex), lty=1)

layout(matrix(1:9, 3), respect=TRUE)
foo = apply(genes, 1, analyse_meth, sunexp_data, sunexp_platform, cols, PLOT=TRUE)

```

