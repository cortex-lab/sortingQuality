

function [clusterIDs, unitQuality, contaminationRate] = maskedClusterQualitySparse(clu, fet, fetInds, fetNchans)
% - clu is 1 x nSpikes
% - fet is nSpikes x nPCsPerChan x nInclChans
% - fetInds is nClusters x nInclChans (sorted in descending order of
% relevance for this template)
% - fetN is an integer, the number of features to

if nargin<4
    fetNchans = min(4, size(fetInds,2)); % number of channels to use
end

fetN = fetNchans*size(fet,2); % now number of features total
% fet = reshape(fet, size(fet,1), []);

N = numel(clu);
assert(fetNchans <= size(fet, 3) && size(fet, 1) == N , 'bad input(s)')

clusterIDs = unique(clu);
unitQuality = zeros(size(clusterIDs));
contaminationRate = zeros(size(clusterIDs));


for c = 1:numel(clusterIDs)
    
    theseSp = clu==clusterIDs(c);
    n = sum(theseSp); % #spikes in this cluster
    if n < fetN || n >= N/2
        % cannot compute mahalanobis distance if less data points than
        % dimensions or if > 50% of all spikes are in this cluster
        unitQuality(c) = 0;
        contaminationRate(c) = NaN;
        continue
    end
    
    fetThisCluster = reshape(fet(theseSp,:,1:fetN), n, []);
    
    % now we need to find other spikes that exist on the same channels
    spikesThisChan = zeros(N,fetNchans,size(fet,2));
    otherFet = fet(~theseSp,:,:);
%     otherFetInds = fetInds(~theseSp,:);
    for chInd = 1:fetNchans
        
        hasThisChan = otherFetInds==fetInds(
        
        spikesThisChan(chInd,:) = clu~=clusterIDs(c) & any(bsxfun(@eq, fetInds
    
    [uQ, cR] = maskedClusterQualityCore(fetThisCluster, fet(i, bestFeatures));
    
    unitQuality(c) = uQ;
    contaminationRate(c) = cR;
    
    fprintf('cluster %d: \t%2.1f\t%1.2f\n', c, unitQuality(c), contaminationRate(c)); % comment to suppress printing out the intermediate results
end
