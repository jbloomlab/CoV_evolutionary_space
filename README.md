# The space of realized and possible coronavirus evolution

Analysis of evolutionary and functional data on SARS-CoV-2 and related coronavirus receptor-binding domains.

## Summary of repository contents

The [input data](./data) in this analysis includes:

 1. Deep mutational scanning measurements of the impact of mutations on ACE2 receptor binding, as measured in [this study](https://www.biorxiv.org/content/10.1101/2022.02.24.481899v1). See [this visualization](https://jbloomlab.github.io/SARS-CoV-2-RBD_DMS_variants/RBD-heatmaps/) for interactive data access and a link to the raw data.

 2. Counts of observations of each mutation on all high-quality, human-derived SARS-CoV-2 sequences present on GISAID as of 27 September, 2021.

 3. An alignment of SARS-related coronavirus RBD sequences, including those from clades capable of ACE2 utilization as described in [this study](https://www.nature.com/articles/s41586-022-04464-z).

The [Rmd script](./analyze_available_space.md) in this repository processes these data and sets thresholds for defining mutations that are tolerated functionally or observed at reasonable frequencies during pandemic SARS-CoV-2 or longer-term sarbecovirus evolution. These metrics are annotated as indicator variables in the output file [here](./results/mutation-specs.csv)
