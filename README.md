# DiffDist
A hierarchical model for RNA-seq data based on the negative binomial distribution, with tests for differential expression, dispersion and distribution. Inference is via an adaptive Markov chain Monte Carlo sampling procedure.

Between-group differences in posterior samples are used to create tail probabilities for differential expression and dispersion based on highest posterior density intervals. Differential expression and dipsersion can be tested for a given minimum log2 fold change.

Differential distribution test is via a mixture model. Inference is based on posterior samples of the relative probability of belonging to each of two mixture components: differential distribution (mean and/or dispersion differs between groups) or no differential distribution (mean and dispersion are the same in both groups). A posterior probability threshold for calling a gene as differentially distributed can be set using a Bayesian false discovery rate control procedure or by using a posterior estimate of the proportion of genes in the sample that are differentially distributed.

Normalisation of RNA-seq counts is not performed as part of this procedure -- normalised counts must be provided. The functions take count matrices with samples in rows and genes in columns -- i.e. the opposite of the conventional format. This will be fixed at some point, but for now count matrices need to be transformed relative to the usual format.
