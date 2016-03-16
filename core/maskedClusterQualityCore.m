
function [unitQuality, contaminationRate] = maskedClusterQualityCore(fetThisCluster, fetOtherClusters)
% fetThisCluster and fetOtherClusters are size [nSpikes, nFeatures]


n = size(fetThisCluster,1);
nOther = size(fetOtherClusters,1);
nFet = size(fetThisCluster,2);

if nOther > n && n>nFet
    % Mahalanobis distance of each of the spikes from present cluster,
      % using only the best fetN dimensions:  
    md = mahal(fetOtherClusters, fetThisCluster);
    md = sort(md);

    mdSelf = mahal(fetThisCluster, fetThisCluster);
    mdSelf = sort(mdSelf);

    unitQuality = md(n);
    contaminationRate = 1-tippingPoint(mdSelf, md)/numel(mdSelf);
else
    unitQuality = 0;
    contaminationRate = NaN;
end

end


function pos = tippingPoint(x,y)
% Input: x, y  are sorted ascending arrays of positive numbers
% Output: minimal pos s.t. sum(x > x(pos)) <= sum(y < x(pos))

% algorithm here is to sort x and y together, and determine the indices of
% x in this sorted list (call this xInds). Then, xInds-(1:length(xInds))
% will be the number of y's that are less than that value of x.

nX = numel(x);
[~, inds] = sort([x;y]);
[~, inds] = sort(inds);
xInds = inds(1:nX);

pos = find(nX:-1:1 < xInds'-(1:nX), 1)-1;

if isempty(pos)
    % not a single "other" spike was nearer the cluster than the furthest
    % in-cluster spike
    pos = nX; % will give contaminationRate = 0;
end

end