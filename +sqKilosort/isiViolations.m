

function isiV = isiViolations(resultsDirectory)

%% Precompute the locationsn of files to be loaded
spikeClustersPath = fullfile(resultsDirectory,'spike_clusters.npy');
spikeTemplatesPath = fullfile(resultsDirectory,'spike_templates.npy');
spikeTimesPath= fullfile(resultsDirectory,'spike_times.npy');
paramsPath= fullfile(resultsDirectory,'params.py');

%% 

refDur = 0.0015;
minISI = 0.0005;

fprintf(1, 'loading data for ISI computation\n');
if exist(spikeClustersPath)
    spike_clusters = readNPY(spikeClustersPath);
else
    spike_clusters = readNPY(spikeTemplatesPath);
end
spike_clusters = spike_clusters + 1; % because in Python indexes start at 0

spike_times = readNPY(spikeTimesPath);
params = readKSparams(paramsPath);
spike_times = double(spike_times)/params.sample_rate;

fprintf(1, 'computing ISI violations\n');

clusterIDs = unique(spike_clusters);
isiV = zeros(1,numel(clusterIDs));
for c = 1:numel(clusterIDs)
    
    [fpRate, numViolations] = ISIViolations(spike_times(spike_clusters==clusterIDs(c)), minISI, refDur);
    isiV(c) = fpRate;
    nSpikes = sum(spike_clusters==clusterIDs(c));    
    fprintf(1, 'cluster %3d: %d viol (%d spikes), %.2f estimated FP rate\n', clusterIDs(c), numViolations, nSpikes, fpRate);
    
end
