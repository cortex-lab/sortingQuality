# sortingQuality
Functions to assess quality of spike sorting results.

The outputs of these functions may serve as a guide to determining how "well isolated" or how contaminated your sorted spikes may be. We do not make any general recommendations about what values returned by these functions may consitute a "good" neuron, and each measure may have its limiations. Ultimately, for any particular experimental result, we recommend making a plot of the value of that result against the sorting quality measure - if the result goes to zero as sorting quality metrics increase, then it likely reflects an artifact of poor spike sorting. 

Ideally, this repository will contain functions in both python and matlab. 

## Function list

### ISIViolations

Computes estimated false positive rate of your spike train, based on the rate of refractory period violations. 

Inputs: spike times, in seconds [1 x nSpikes]; minimum possible inter-spike interval given your spike detection algorithm, in seconds; refractory period duration, in seconds;

Output: estimated rate of spikes that come from another neuron besides the primary one. 

Notes: This analysis assumes the primary neuron and any other neurons are completely uncorrelated in their spike times. It will be highly unreliable in cases when this assumption is strongly violated, e.g. when two visual cortical neurons have different orientation preference and, when driven with oriented stimuli, rarely spike together due to the stimulus set. This analysis will be unreliable for neurons with low spike counts; in particular, if the neuron has a low spike rate, then zero refractory violations may be a statistically likely outcome even in the face of strong contamination. 
