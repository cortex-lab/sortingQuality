

function [fpRate, numViolations] = ISIViolations(spikeTrain, minISI, refDur)
% computes an estimated false positive rate of the spikes in the given
% spike train. You probably don't want this to be greater than a few
% percent of the spikes. 
%
% - minISI [sec] = minimum *possible* ISI (based on your cutoff window); likely 1ms
% or so. It's important that you have this number right.
% - refDur [sec] = estimate of the duration of the refractory period; suggested value = 0.002.
% It's also important that you have this number right, but there's no way
% to know... also, if it's a fast spiking cell, you could use a smaller
% value than for regular spiking.

totalRate = length(spikeTrain)/spikeTrain(end);
numViolations = sum(diff(spikeTrain) <= refDur);

% Spikes occurring under the minimum ISI are by definition false positives;
% their total rate relative to the rate of the neuron overall gives the
% estimated false positive rate.
% NOTE: This does not use the equation from Dan Hill's paper (see below)
% but instead uses the equation from his UltraMegaSort
violationTime = 2*length(spikeTrain)*(refDur-minISI); % total time available for violations - 
                                                    % there is an available
                                                    % window of length
                                                    % (refDur-minISI) after
                                                    % each spike.
violationRate = numViolations/violationTime;
fpRate = violationRate/totalRate;

if fpRate>1
    % it is nonsense to have a rate >1, however having such large rates
    % does tell you something interesting, namely that the assumptions of
    % this analysis are failing!
    fpRate = 1; 
end


% Here is with Dan Hill's equation (J. Neurosci, 2012)
%  NOTE: actually, this method will return imaginary results for large-ish
%  numbers of violations. 
% tauR = refDur;
% tauC = minISI;
% T = spikeTrain(end);
% N = length(spikeTrain);
% R = sum(diff(spikeTrain) <= tauR);
% 
% k = 2*(tauR-tauC)*N^2;
% rts = roots([k -k R*T]);
% 
% fpRate = min(rts);
% numViolations = R;
end