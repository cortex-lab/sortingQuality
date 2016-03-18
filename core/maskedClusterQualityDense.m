


function [clusterIDs, unitQuality, contaminationRate] = maskedClusterQualityDense(clu, fet, fmask, fetN)
% function [clusterIDs, unitQuality, contaminationRate] = maskedClusterQualityDense(clu, fet, fmask, fetN)
% Inputs:
% clu   N-by-1 vector which assigns each spike to a cluster (here N is total no. of spikes)
% fet   N-by-(total #features) matrix of features of each spike
% fmask N-by-(total #features) matrix of masks for each feature of every spike
% fetN  (optional, 12 by default). # of features to take into consideration for quality computation
% Output: clusterIDs, unitQuality, contaminationRate
% Note that for a cluster that contains over 50% of all spikes, or less than fetN spikes, quality is not defined (0)
%
% Code by M. Okun, some stylistic changes by N. Steinmetz

if nargin < 4
  fetN = 12;
end
N = numel(clu);
assert(size(fet, 2) == size(fmask, 2) && fetN <= size(fet, 2) && ...
  size(fet, 1) == N && size(fmask, 1) == N, 'bad input(s)')

clusterIDs = unique(clu);
unitQuality = zeros(size(clusterIDs));
contaminationRate = zeros(size(clusterIDs));

for c = 1:numel(clusterIDs)
  n = sum(clu == clusterIDs(c)); % #spikes in this cluster
  if n < fetN || n >= N/2
    % cannot compute mahalanobis distance if less data points than
    % dimensions or if > 50% of all spikes are in this cluster
    unitQuality(c) = 0;
    contaminationRate(c) = NaN;
    continue
  end
  
  [~, bestFeatures] = sort(mean(fmask(clu == clusterIDs(c), :)), 'descend');
  bestFeatures = bestFeatures(1:fetN); % the best fetN features for this cluster
    
  % We don't want to take into consideration spikes that have absolutely
  % nothing to do with these dimensions:
  i = (clu ~= clusterIDs(c)) & sum(fmask(:, bestFeatures), 2) > 0;
    
  [uQ, cR] = maskedClusterQualityCore(fet(clu == clusterIDs(c), bestFeatures), fet(i, bestFeatures));
  
  unitQuality(c) = uQ;
  contaminationRate(c) = cR;
  
  fprintf('cluster %d: \t%2.1f\t%1.2f\n', clusterIDs(c), unitQuality(c), contaminationRate(c)); % comment to suppress printing out the intermediate results
end

