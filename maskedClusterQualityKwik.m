
function [clusterIDs, unitQuality, contaminationRate] = maskedClusterQualityKwik(baseFilename, shankIndex)
% if ispc
%   baseFilename  = '\\basket.cortexlab.net\M150218_NS1LAV\20150601\20150601_all';
% else
%   baseFilename  = '/data/nick/M150218_NS1LAV/20150601/20150601_all';
% end

clusterNames =  hdf5read([baseFilename '.kwik'], ['/channel_groups/' num2str(shankIndex) '/spikes/clusters/main']);
featuresAndMasks = hdf5read([baseFilename '.kwx'],  ['/channel_groups/' num2str(shankIndex) '/features_masks']); %  Size is [2, nChan*3, nSpikes]

clusterNames = double(clusterNames);
features = squeeze(featuresAndMasks(1,:,:))';
masks = squeeze(featuresAndMasks(2,:,:))';
clear featuresAndMasks

[clusterIDs, unitQuality, contaminationRate] = maskedClusterQualityDense(clusterNames, features, masks);