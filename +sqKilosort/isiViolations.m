

function isiV = isiViolations(resultsDirectory)

refDur = 0.002;
minISI = 0.0005;

fprintf(1, 'loading data for ISI computation\n');
if exist([resultsDirectory 'spike_clusters.npy'])
    spike_clusters = readNPY([resultsDirectory 'spike_clusters.npy']);
else
    spike_clusters = readNPY([resultsDirectory 'spike_templates.npy']);
end
spike_clusters = spike_clusters + 1; % because in Python indexes start at 0

spike_times = readNPY([resultsDirectory 'spike_times.npy']);
params = readKSparams([resultsDirectory 'params.py']);
spike_times = double(spike_times)/params.sample_rate;

fprintf(1, 'computing ISI violations\n');

clusterIDs = unique(spike_clusters);
isiV = zeros(1,numel(clusterIDs));
for c = 1:numel(clusterIDs)
    
    [fpRate, numViolations] = ISIViolations(spike_times(spike_clusters==clusterIDs(c)), minISI, refDur);
    isiV(c) = fpRate;
    nSpikes = sum(spike_clusters==clusterIDs(c));    
    fprintf(1, 'cluster %d: %d viol (%d spikes), %.2f estimated FP rate\n', clusterIDs(c), numViolations, nSpikes, fpRate);
    
end
