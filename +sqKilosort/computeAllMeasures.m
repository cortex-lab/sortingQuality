

function [cgs, uQ, cR, isiV] = computeAllMeasures(resultsDirectory)

if exist(fullfile(resultsDirectory, 'cluster_groups.csv'))
    [cids, cgs] = readClusterGroupsCSV(fullfile(resultsDirectory, 'cluster_groups.csv'));
elseif exist(fullfile(resultsDirectory, 'spike_clusters.npy'))
    clu = readNPY(fullfile(resultsDirectory, 'spike_clusters.npy'));
    cgs = 3*ones(size(unique(clu))); % all unsorted
else
    clu = readNPY(fullfile(resultsDirectory, 'spike_templates.npy'));
    cgs = 3*ones(size(unique(clu))); % all unsorted
end

[cids, uQ, cR] = sqKilosort.maskedClusterQuality(resultsDirectory);

isiV = sqKilosort.isiViolations(resultsDirectory);

