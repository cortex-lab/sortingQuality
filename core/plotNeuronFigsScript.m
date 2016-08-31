ksRoot = 'J:\M160731_MOEC\2016-08-28';
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
d = dir(params.filename); nSamp = d.bytes/2/sp.n_channels_dat;
params.dataSize = [sp.n_channels_dat nSamp];
params.chanMap = readNPY('channel_map.npy');
params.Fs = sp.sample_rate;
params.xcoords = sp.xcoords; params.ycoords = sp.ycoords;
params.plotDistance = 100;
params.nWFsToLoad = 1000;
params.nWFsToPlot = 100;
params.gain = 0.6/512/500*1e6; % raw file units to uV
params.nPCsToPlot = 50000;
sparsePCfeat = sparsePCs(pcFeat, pcFeatInd, sp.spikeTemplates);
for q = 1%:length(inclCID)
    thisClu = inclCID(q);

    [v, i] = countUnique(sp.spikeTemplates(clu==thisClu));
    thisTempID = v(find(i==max(i),1));
    template = squeeze(sp.temps(thisTempID+1,:,:));

    figHand = neuronFig(thisClu, st, clu, sparsePCfeat, [], params);
    set(figHand, 'Position', [-1890         -59        1810        1031]);
%     saveFig(figHand, fullfile(figDir, sprintf('/cluster%d', thisClu)))
%     close(figHand); clear figHand
%     makeClusterFigsWebsite(figDir)
end