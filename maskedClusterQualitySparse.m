

function [clusterIDs, unitQuality, contaminationRate] = maskedClusterQualitySparse(clu, fet, fetInds, fetNchans)
% - clu is 1 x nSpikes
% - fet is nSpikes x nPCsPerChan x nInclChans
% - fetInds is nClusters x nInclChans (sorted in descending order of
% relevance for this template)
% - fetN is an integer, the number of features to

if nargin<4
    fetNchans = min(4, size(fetInds,2)); % number of channels to use
end
nFetPerChan = size(fet,2);
fetN = fetNchans*nFetPerChan; % now number of features total
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
    
    fetThisCluster = reshape(fet(theseSp,:,1:fetNchans), n, []);
    
    % now we need to find other spikes that exist on the same channels
    theseChans = fetInds(c,1:fetNchans);
    for f = 1:fetNchans
        thisChanInds = fetInds==theseChans(f);
        [chanInds,clustWithThisChan] = find(thisChanInds');
        chanInds = chanInds(clustWithThisChan~=c);
        clustWithThisChan = clustWithThisChan(clustWithThisChan~=c);                
        
        otherSpikes = ismember(clu, clusterIDs(clustWithThisChan));
        nOtherSpikes = sum(otherSpikes);
        
        fetOtherClusters = zeros(nOtherSpikes, nFetPerChan, fetNchans);
        nInd = 1;              
        
        for t = 1:numel(clustWithThisChan)            
            thisCfetInd = chanInds(t);
            theseOtherSpikes = clu==clusterIDs(clustWithThisChan(t));
            fetOtherClusters(nInd:nInd+sum(theseOtherSpikes)-1,:,f) = ...
                fet(theseOtherSpikes,:,thisCfetInd);
            nInd = nInd+sum(theseOtherSpikes);
        end
        
    end
    
    fetOtherClusters = reshape(fetOtherClusters, size(fetOtherClusters,1), []);
    
    [uQ, cR] = maskedClusterQualityCore(fetThisCluster, fetOtherClusters);
    
    unitQuality(c) = uQ;
    contaminationRate(c) = cR;
    
    fprintf('cluster %d: \t%2.1f\t%1.2f\n', c, unitQuality(c), contaminationRate(c)); % comment to suppress printing out the intermediate results
end
