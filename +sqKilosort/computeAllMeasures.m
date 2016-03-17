

function [cgs, uQ, cR, isiV] = computeAllMeasures(resultsDirectory)

if exist(fullfile(resultsDirectory, 'cluster_groups.csv'))
    [cids, cgs] = readClusterGroupsCSV(fullfile(resultsDirectory, 'cluster_groups.csv'));
else
    clu = readNPY(fullfile(resultsDirectory, 'spike_clusters.npy'));
    cgs = 3*ones(size(unique(clu))); % all unsorted
end

[cids, uQ, cR] = sqKilosort.maskedClusterQuality(resultsDirectory);

isiV = sqKilosort.isiViolations(resultsDirectory);

