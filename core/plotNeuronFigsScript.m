% ksRoot = 'G:\data\Hopkins\2016-07-22';
% ksRoot = 'J:\M160731_MOEC\2016-08-28'
% ksRoot = 'J:\Hopkins\20160722'
% ksRoot = 'J:\Eijkman\2016-05-21\M2';
% ksRoot = 'J:\Eijkman\2016-05-21\visctx';
% ksRoot = 'J:\Loewi\2016-04-20\frontal';
% ksRoot = 'J:\SS060\20160520\SC';
% ksRoot = 'J:\SS061\20160511\';
% ksRoot = 'J:\Loewi\2016-04-20\posterior';
ksRoot = 'J:\M160731_MOEC\2016-09-16';

loadPars.loadPCs = true;
sp = loadKSdir(ksRoot, loadPars);

inclCID = sp.cids(sp.cgs==2);
st = sp.st;
clu = sp.clu;

pcFeat = sp.pcFeat;
pcFeatInd = sp.pcFeatInd;
spikeAmps = sp.tempScalingAmps;
figDir = fullfile(ksRoot, 'figs'); 
if ~exist(figDir,'dir'); mkdir(figDir); end

params.dataType = sp.dtype;
params.filename = sp.dat_path;
d = dir(fullfile(ksRoot,params.filename)); 

% params.filename = 'V:\nick\Eijkman\2016-05-21\Eijkman_20160521_visctx_g0_t0.imec_AP_CAR.bin';
% d = dir(params.filename); 

nSamp = d.bytes/2/sp.n_channels_dat;
params.dataSize = [sp.n_channels_dat nSamp];
params.chanMap = readNPY(fullfile(ksRoot, 'channel_map.npy'));
params.Fs = sp.sample_rate;
params.xcoords = sp.xcoords; params.ycoords = sp.ycoords;
params.plotDistance = 100;
params.nWFsToLoad = 1000;
params.nWFsToPlot = 100;
params.gain = 0.6/512/500*1e6; % raw file units to uV
params.nPCsToPlot = 50000;
% params.highlightRegions = [132 1104];

%% extract median WFs (just once)
cd(ksRoot)
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
theseCR = cR(cgs==2);
theseID = uQ(cgs==2);

for q = 29:length(inclCID)
    clusterID = inclCID(q);

    stats.medWF = squeeze(medWFs(inclCID==clusterID,:,:))';
    stats.isiContamination = theseISI(inclCID==clusterID);
    stats.isoDistance = theseID(inclCID==clusterID);
    stats.mahalContamination = theseCR(inclCID==clusterID);
%     stats.isiContamination = 0;
%     stats.isoDistance = 0;
    figHand = neuronFig(clusterID, st, clu, sparsePCfeat, spikeAmps, stats, params);
    set(figHand, 'Position', [-1890         -59        1810        1031]);
    saveFig(figHand, fullfile(figDir, sprintf('/cluster%d', clusterID)))
    close(figHand); clear figHand
    
end

%%

[~, spikeDepths, ~, ~, ~, ~, ~] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
cluDepths = clusterAverage(sp.clu(ismember(sp.clu,sp.cids(sp.cgs==2))), spikeDepths(ismember(sp.clu,sp.cids(sp.cgs==2))));
makeClusterFigsWebsite(figDir, 'sortByDepth.html', cluDepths);

cluAmps = max(max(medWFs,[],3)-min(medWFs,[],3), [], 2);
makeClusterFigsWebsite(figDir, 'sortByAmp.html', cluAmps);

[~, FRs] = countUnique(sp.clu(ismember(sp.clu, sp.cids(sp.cgs==2))));
makeClusterFigsWebsite(figDir, 'sortByFR.html', FRs);

makeClusterFigsWebsite(figDir, 'sortByIsoDist.html', uQ(cgs==2));
