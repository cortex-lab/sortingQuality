# sortingQuality
Functions to assess quality of spike sorting results.

The outputs of these functions may serve as a guide to determining how "well isolated" or how contaminated your sorted spikes may be. We do not make any general recommendations about what values returned by these functions consitute a "good" neuron, and each measure has its limiations. Ultimately, for any particular experimental result, we recommend making a plot of the value of that result against any sorting quality measures - if the result disappears as sorting quality metrics improve, then it likely reflects an artifact of poor spike sorting. 

Ideally, this repository will contain functions in both python and matlab. 

## Function list

### ISIViolations

Computes estimated false positive rate of your spike train, based on the rate of refractory period violations. Method from Hill et al., J. Neurosci, 2011. See [their function rpv_contamination](https://github.com/rheitz1/Mat_Code/blob/master/UltraMegaSort/quality_measures/rpv_contamination.m). 

Inputs: spike times, in seconds [1 x nSpikes]; minimum possible inter-spike interval given your spike detection algorithm, in seconds; refractory period duration, in seconds;

Output: estimated rate of spikes that come from another neuron besides the primary one. 

Notes: This analysis assumes the primary neuron and any other neurons are completely uncorrelated in their spike times. It will be highly unreliable in cases when this assumption is strongly violated, e.g. when two visual cortical neurons have different orientation preference and, when driven with oriented stimuli, rarely spike together due to the stimulus set. This analysis will be unreliable for neurons with low spike counts; in particular, if the neuron has a low spike rate, then zero refractory violations may be a statistically likely outcome even in the face of strong contamination. 

### maskedClusterQuality

Computes "isolation distance". If there are more than four channels on the probe, uses a modified calculation with only the top four channels per cluster, to avoid dimensionality problems. In addition to returning Isolation Distance, this method also returns a "contamination rate" which is the proportion of spikes inside the cluster boundary that aren't from the cluster (false positive rate), if you set the cluster boundary at a mahalanobis distance such that there are equal false positives and false negatives. Method adapted from Schmitzer-Torber et al., Neuroscience, 2005. 

### GaussianContamination

Estimates false positive and false negative rates for each cluster by fitting N-dim gaussians, sampling from these distributions, classifying the sampled points, and computing the error from these samples. Method from Tolias et al., J. Neurophys, 2007. 

### MultiClustering

Computes estimate of quality by re-running clustering multiple times and estimating consistency between runs.

## References

Hill, D. N., Mehta, S. B., & Kleinfeld, D. (2011). Quality metrics to accompany spike sorting of extracellular signals. Journal of Neuroscience, 31(24), 8699–705. doi:10.1523/JNEUROSCI.0971-11.2011

Schmitzer-Torbert, N., Jackson, J., Henze, D. A., Harris, K. D., & Redish, a. D. (2005). Quantitative measures of cluster quality for use in extracellular recordings. Neuroscience, 131(1), 1–11. doi:10.1016/j.neuroscience.2004.09.066

Tolias, A. S., Ecker, A. S., Siapas, A. G., Hoenselaar, A., Keliris, G. A., & Logothetis, N. K. (2007). Recording chronically from the same neurons in awake, behaving primates. Journal of Neurophysiology, 98(6), 3780. doi:10.1152/jn.00260.2007
