


function isiV = isiViolations(baseFilename, shankIndex)

refDur = 0.002;
minISI = 0.0005;

kwikFilename = [baseFilename '.kwik'];

fprintf(1, 'loading data for ISI computation\n');
spike_clusters = h5read(kwikFilename, ['/channel_groups/' num2str(shankIndex) '/spikes/clusters/main']);

spikeSamps = h5read(kwikFilename, ['/channel_groups/' num2str(shankIndex) '/spikes/time_samples']);
Fs = double(h5readatt(kwikFilename, '/recordings/0', 'sample_rate'));
spike_times = double(spikeSamps)/Fs;

fprintf(1, 'computing ISI violations\n');

clusterIDs = unique(spike_clusters);
isiV = zeros(1,numel(clusterIDs));
for c = 1:numel(clusterIDs)
    
    [fpRate, numViolations] = ISIViolations(spike_times(spike_clusters==clusterIDs(c)), minISI, refDur);
    isiV(c) = fpRate;
    nSpikes = sum(spike_clusters==clusterIDs(c));    
    fprintf(1, 'cluster %d: %d viol (%d spikes), %.2f estimated FP rate\n', clusterIDs(c), numViolations, nSpikes, fpRate);
    
end
