

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

fprintf('%12s\tQuality\tContamination\n', 'ID'); % comment to suppress printing out the intermediate results
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
%     for f = 1:fetNchans
%         thisChanInds = fetInds==theseChans(f);
%         [chanInds,clustWithThisChan] = find(thisChanInds');
%         chanInds = chanInds(clustWithThisChan~=c);
%         clustWithThisChan = clustWithThisChan(clustWithThisChan~=c);                
%         
%         otherSpikes = ismember(clu, clusterIDs(clustWithThisChan));
%         nOtherSpikes = sum(otherSpikes);
%         
%         fetOtherClusters = zeros(nOtherSpikes, nFetPerChan, fetNchans);
%         nInd = 1;              
%         
%         for t = 1:numel(clustWithThisChan)            
%             thisCfetInd = chanInds(t);
%             theseOtherSpikes = clu==clusterIDs(clustWithThisChan(t));
%             fetOtherClusters(nInd:nInd+sum(theseOtherSpikes)-1,:,f) = ...
%                 fet(theseOtherSpikes,:,thisCfetInd);
%             nInd = nInd+sum(theseOtherSpikes);
%         end
%         
%     end

    % for each other cluster, determine whether it has at least one of
    % those channels. If so, add its spikes, with its features put into the
    % correct places
    nInd = 1; fetOtherClusters = zeros(0,size(fet,2),fetNchans);
    for c2 = 1:numel(clusterIDs)
        if c2~=c
            chansC2Has = fetInds(c2,:);
            for f = 1:length(theseChans)
                
                if ismember(theseChans(f), chansC2Has)
                    
                    theseOtherSpikes = clu==clusterIDs(c2);
                    thisCfetInd = find(chansC2Has==theseChans(f),1);
                    fetOtherClusters(nInd:nInd+sum(theseOtherSpikes)-1,:,f) = ...
                        fet(theseOtherSpikes,:,thisCfetInd);
                end                
                
            end
            if any(ismember(chansC2Has, theseChans))
                nInd = nInd+sum(theseOtherSpikes);
            end
        end
    end
    
                
                
                
    
    fetOtherClusters = reshape(fetOtherClusters, size(fetOtherClusters,1), []);
    
    [uQ, cR] = maskedClusterQualityCore(fetThisCluster, fetOtherClusters);
    
    unitQuality(c) = uQ;
    contaminationRate(c) = cR;
    
    fprintf('cluster %3d: \t%6.1f\t%6.2f\n', clusterIDs(c), unitQuality(c), contaminationRate(c)); % comment to suppress printing out the intermediate results
    
    if uQ>1000
        keyboard;
    end
    
end
