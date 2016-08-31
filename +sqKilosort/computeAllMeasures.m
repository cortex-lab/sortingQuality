

function [cgs, uQ, cR, isiV] = computeAllMeasures(resultsDirectory)

clusterPath = fullfile(resultsDirectory, 'cluster_groups.csv');
spikeClustersPath = fullfile(resultsDirectory,'spike_clusters.npy');
spikeTemplatesPath = fullfile(resultsDirectory,'spike_templates.npy');

if exist(clusterPath, 'file')
    [cids, cgs] = readClusterGroupsCSV(clusterPath);
elseif exist(spikeClustersPath, 'file')
    clu = readNPY(spikeClustersPath);
    cgs = 3*ones(size(unique(clu))); % all unsorted
else
    clu = readNPY(spikeTemplatesPath);
    cgs = 3*ones(size(unique(clu))); % all unsorted
end

[cids, uQ, cR] = sqKilosort.maskedClusterQuality(resultsDirectory);

isiV = sqKilosort.isiViolations(resultsDirectory);

