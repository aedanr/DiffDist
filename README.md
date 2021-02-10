# DiffDist
https://zenodo.org/badge/latestdoi/337280917

A hierarchical model for RNA-seq data based on the negative binomial distribution, with tests for differential expression, dispersion and distribution. Inference is via an adaptive Markov chain Monte Carlo sampling procedure.

Between-group differences in posterior samples are used to create tail probabilities for differential expression and dispersion based on highest posterior density intervals. Differential expression and dipsersion can be tested for a given minimum log2 fold change.

Differential distribution test is via a mixture model. Inference is based on posterior samples of the relative probability of belonging to each of two mixture components: differential distribution (mean and/or dispersion differs between groups) or no differential distribution (mean and dispersion are the same in both groups). A posterior probability threshold for calling a gene as differentially distributed can be set using a Bayesian false discovery rate control procedure or by using a posterior estimate of the proportion of genes in the sample that are differentially distributed.

Normalisation of RNA-seq counts is not performed as part of this procedure -- normalised counts must be provided. The functions take count matrices with samples in rows and genes in columns -- i.e. the opposite of the conventional format. This will be fixed at some point, but for now count matrices need to be transformed relative to the usual format.

`lognormal_hmm_adaptive_proposals_three_chains_function` contains the main function to run the MCMC algorithm. This function calls the separate functions contained in `lognormal_hmm_three_chains_function` and `lognormal_hmm_one_chain_function`, which rely on the functions in `conditional_posterior_functions_lognormal_hmm`.

`methods_comparisons` contains code used to run the DiffDist tests, along with code used to run comparative assessments against other differential expression, dispersion and distribution methods. The code for the hierarchical model tests here also serves as an example of how to use the output of the MCMC algorithm to perform each test.
