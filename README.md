Each script in this repository corresponds to the following described analyses that are included in the manuscript.

anc_char_est.R: Calculates the tip rates for each of the species included in our analyses.

carn.chrom.analysis.R: chromePlus analysis of the carnivore chromosome number and range size data

carn.discrete.R: find value for discretzation of range size based on our carnivore data

carn.sim.binary.analysis.R: analysis that simulates neutrally evolving traits under chromePlus model to estimate a false positive rate

cor.est.R: compares published range size estimates to those that we estimated for our analyses

data_munging_rangesize.R: scripot that matches up tree data and chromosome number data, as well as discretizes our binary trait

functions.R: functions needed to perform the scalar analysis

getQ.function.R: a function that creates the Q matrix needed before estimating tip rates

process.data.R: processes intial data gathering to keep only those species that we have a tip for in our phylogeny

result_munging_rangesize.R: munges tha results from each of the chromePlus replicates to gather the final posterior dataset 

scalar.analysis.R: performs the scalar analysis for each of the families in our dataset
