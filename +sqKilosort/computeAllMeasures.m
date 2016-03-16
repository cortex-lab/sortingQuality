

function [cgs, uQ, cR, isiV] = computeAllMeasures(resultsDirectory)

[cids, cgs] = sqKilosort.readClusterGroupsCSV(fullfile(resultsDirectory, 'cluster_groups.csv'));

[