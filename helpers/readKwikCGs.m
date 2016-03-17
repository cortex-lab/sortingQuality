

function cgs = readKwikCGs(baseFilename, shankIndex)

kwikFilename = [baseFilename '.kwik'];

spikeClusterIDs = h5read(kwikFilename, ['/channel_groups/' num2str(shankIndex) '/spikes/clusters/main']);

clusterIDs = unique(spikeClusterIDs);
numClusters = numel(clusterIDs);

cgs = zeros(1, numClusters);

for c = 1:numClusters
        
    cgs(c) = h5readatt(kwikFilename, ['/channel_groups/' num2str(shankIndex) '/clusters/main/' num2str(clusterIDs(c)) '/'], 'cluster_group');
    
end

end
