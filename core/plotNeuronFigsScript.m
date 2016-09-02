% ksRoot = 'G:\data\Hopkins\2016-07-22';
% ksRoot = 'J:\M160731_MOEC\2016-08-28'
ksRoot = 'J:\Hopkins\20160722'

loadPars.loadPCs = true;
sp = loadKSdir(ksRoot, loadPars);

inclCID = sp.cids(sp.cgs==2);
st = sp.st;
clu = sp.clu;
pcFeat = sp.pcFeat;
pcFeatInd = sp.pcFeatInd;

figDir = fullfile(ksRoot, 'figs'); 
if ~exist(figDir,'dir'); mkdir(figDir); end

params.dataType = sp.dtype;
params.filename = sp.dat_path;
d = dir(fullfile(ksRoot,params.filename)); nSamp = d.bytes/2/sp.n_channels_dat;
params.dataSize = [sp.n_channels_dat nSamp];
params.chanMap = readNPY(fullfile(ksRoot, 'channel_map.npy'));
params.Fs = sp.sample_rate;
params.xcoords = sp.xcoords; params.ycoords = sp.ycoords;
params.plotDistance = 100;
params.nWFsToLoad = 1000;
params.nWFsToPlot = 100;
params.gain = 0.6/512/500*1e6; % raw file units to uV
params.nPCsToPlot = 50000;

%% extract median WFs (just once)

inclSP = ismember(clu, sp.cids(sp.cgs==2));
medWFs = extractMedianWFs(clu(inclSP), st(inclSP), params.Fs, params.filename, ...
    params.dataType, params.dataSize, params.chanMap, params.gain);

save(fullfile(ksRoot, 'medWFs.mat'), 'medWFs');

%% otherwise load them

load(fullfile(ksRoot, 'medWFs.mat'))

%% compute cluster quality stats (just once)

[cgs, uQ, cR, isiV] = sqKilosort.computeAllMeasures(ksRoot);

save(fullfile(ksRoot, 'clusterQualityMetrics.mat'), 'cgs', 'uQ', 'cR', 'isiV');

%%

sparsePCfeat = sparsePCs(pcFeat, pcFeatInd, sp.spikeTemplates, 2, 10);

%%
figDir = 'V:\www\data.cortexlab.net\singlePhase3\figs';

theseISI = isiV(cgs==2);
theseID = uQ(cgs==2);

for q = 2:length(inclCID)
    clusterID = inclCID(q);

    stats.medWF = squeeze(medWFs(inclCID==clusterID,:,:))';
    stats.isiContamination = theseISI(inclCID==clusterID);
    stats.isoDistance = theseID(inclCID==clusterID);
%     stats.isiContamination = 0;
%     stats.isoDistance = 0;
    figHand = neuronFig(clusterID, st, clu, sparsePCfeat, stats, params);
    set(figHand, 'Position', [-1890         -59        1810        1031]);
    saveFig(figHand, fullfile(figDir, sprintf('/cluster%d', clusterID)))
    close(figHand); clear figHand
    makeClusterFigsWebsite(figDir)
end