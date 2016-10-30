

function f = neuronFig(clusterID, spikeTimes, clu, sparsePCfeat, spikeAmps, stats, params)
% Makes a plot of relevant stats/figures for a neuron
%
%
% stats has:
% - isoDistance
% - isiContamination
%
% params has: 
% - chanMap
% - xcoords, ycoords
% - "dataType" and "dataSize" of raw file
% - raw "filename"
% - Fs


if isempty(stats)
    stats.isiContamination = NaN;
    stats.isoDistance = NaN;
end

%%

Fs = params.Fs;
wfWin = -round(0.5*Fs/1000):round(1.25*Fs/1000); nWFsamps = numel(wfWin);
nChInFile = params.dataSize(1);
nSamp = params.dataSize(2);
nCh = numel(params.chanMap);
xc = params.xcoords; yc = params.ycoords;

%%

theseST = spikeTimes(clu==clusterID);

% Extract waveforms from this neuron and from background
params.nWFsToLoad = params.nWFsToPlot; % computing meanWF happens elsewhere now, so only load what we'll plot
nWFsToLoad = min(params.nWFsToLoad, length(theseST));
extractST = round(theseST(randperm(length(theseST), nWFsToLoad))*Fs);


extractBckg = randperm(nSamp-numel(wfWin)-2, params.nWFsToPlot)-wfWin(1); % *samples* of background spikes

mmf = memmapfile(params.filename, 'Format', {params.dataType, [nChInFile nSamp], 'x'});

theseWF = zeros(nWFsToLoad, nCh, nWFsamps);
bckgWF = zeros(params.nWFsToPlot, nCh, nWFsamps);
for i=1:nWFsToLoad
     tempWF = mmf.Data.x(1:nChInFile,extractST(i)+wfWin(1):extractST(i)+wfWin(end));
     theseWF(i,:,:) = tempWF(params.chanMap+1,:);
end

for i = 1:params.nWFsToPlot %
     tempWF = mmf.Data.x(1:nChInFile,extractBckg(i)+wfWin(1):extractBckg(i)+wfWin(end));
     bckgWF(i,:,:) = tempWF(params.chanMap+1,:);
     
end

% medWF = squeeze(median(double(theseWF),1))';
% medWFuV = medWF.*params.gain;
medWFuV = stats.medWF;

%%

% choose channels to plot

chanAmps = max(medWFuV)-min(medWFuV);
maxChan = find(chanAmps==max(chanAmps),1);
maxXC = xc(maxChan); maxYC = yc(maxChan);
chanDistances = ((xc-maxXC).^2 + (yc-maxYC).^2).^0.5;
chansToPlot = chanDistances<params.plotDistance;

wfAmp = max(chanAmps);

bckgMaxChan = double(squeeze(bckgWF(:,maxChan,:)));
% snr = wfAmp./std(bckgMaxChan(:));
snr = wfAmp./median(abs(bckgMaxChan(:))/0.6745); % RQQ method

%% compute PC stuff

thesePCs = sparsePCfeat(clu==clusterID,:);
meanPC = mean(thesePCs);
[~, ii] = sort(abs(meanPC), 2, 'descend');
topChans = ii(1:2);

% Method 1: just pick the top two channels for this cluster
otherSpikesIncl = sparsePCfeat(:,topChans(1))~=0 & ...
    sparsePCfeat(:,topChans(2))~=0 & clu~=clusterID;
otherSpikesPCs = sparsePCfeat(otherSpikesIncl, topChans);
np = min(size(otherSpikesPCs,1), params.nPCsToPlot);
otherPCsToPlotInds = randperm(size(otherSpikesPCs,1), np);
otherPCsToPlot = otherSpikesPCs(otherPCsToPlotInds,:);
thesePCsToPlot = thesePCs(:,topChans);

% % Method 2: figure out the top two PCs *of the PCs* of this neuron, project
% % all other spikes onto those
% pcInclChans = full(meanPC)~=0;
% otherSpikesIncl = sum(sparsePCfeat(:,pcInclChans)~=0,2)>0; % those spikes with non-zero values on at least one included channel
% otherSpikesPCs = sparsePCfeat(otherSpikesIncl, pcInclChans);
% 
% thesePCsIncl = thesePCs(:,pcInclChans);
% [coeff, score, latent, tsquared, explained, mu] = pca(full(thesePCsIncl), 'Centered', false);
% thesePCsToPlot = score(:,1:2);
% otherPCsToPlotInds = randperm(size(otherSpikesPCs,1), params.nPCsToPlot);
% % project the other spikes to be plotted onto these new vectors
% otherSpikesToPlotPCs = otherSpikesPCs(otherPCsToPlotInds,:);
% otherPCsToPlot = otherSpikesToPlotPCs*coeff(:,1:2);


% % Method 3: one dimension is the mean PC of the cluster. Other dimension is
% % the direction between this cluster and the nearest cluster. 
% pcInclChans = full(meanPC)~=0;
% otherSpikesIncl = sum(sparsePCfeat(:,pcInclChans)~=0,2)>0 & clu~=clusterID; % those spikes with non-zero values on at least one included channel
% otherSpikesClu = clu(otherSpikesIncl);
% otherSpikesPCs = sparsePCfeat(otherSpikesIncl, pcInclChans); % nOtherSpikes x nIncludedPC
% thesePCsIncl = thesePCs(:,pcInclChans); % nSpikes x nIncludedPC
% meanPCincl = meanPC(pcInclChans); % 1 x nIncludedPC
% 
% theseScore1 = thesePCsIncl*meanPCincl'; % nSpikes x 1
% otherScore1 = otherSpikesPCs*meanPCincl'; % nOtherSpikes x 1
% 
% % this would be the linear discriminant if incorporated covariance
% % matrices. But not sure that makes sense since there are so many other
% % clusters...?
% % direction2 = mean(theseResidPC)-mean(otherResidPC); % 1 x nIncludedPC 
% % direction2 = mean(thesePCsIncl)-mean(otherSpikesPCs); % 1 x nIncludedPC
% 
% otherCID = unique(otherSpikesClu);
% nOtherClu = length(otherCID);
% otherMeanPCs = zeros(nOtherClu, sum(pcInclChans));
% for q = 1:sum(pcInclChans)
%     otherMeanPCs(:,q) = clusterAverage(otherSpikesClu, otherSpikesPCs(:,q));
% end
% similarity = meanPCincl*otherMeanPCs';
% closestInd = find(similarity==max(similarity),1);
% direction2 = meanPCincl-otherMeanPCs(closestInd,:);
% 
% theseScore2 = thesePCsIncl*direction2'; % nSpikes x 1
% otherScore2 = otherSpikesPCs*direction2'; % nOtherSpikes x 1
% 
% thesePCsToPlot = [theseScore1 theseScore2];
% otherPCsToPlot = [otherScore1 otherScore2];

%%
f = figure;
set(f, 'Color', 'w');

% neuronColor = [0.8500    0.3250    0.0980]; % the orange you get with second element in default color order
neuronColor = [ 0    0.4470    0.7410]; % blue that's first default
otherColor = [0.4660    0.6740    0.1880];

%% plot of location on probe
subplot(4,5,[1 6 11 16])
plotAsProbe(-chanAmps, xc, yc, colormap(gray), 16, 40)
hold on;
minx = min(xc)-8; maxx = max(xc)+8; miny = min(yc)-20; maxy = max(yc)+20;
plot([minx maxx maxx minx minx], [maxy maxy miny miny maxy], 'k');
axis off
title('neuron position')
% xlabel('space (µm)')
% ylabel('space (µm)')

%% plot of waveforms, zoomed around location
subplot(4,5, [2 3 7 8 12 13]);
xPlot = xc(chansToPlot);
yPlot = yc(chansToPlot);
p.LineWidth = 0.01; p.alpha = 0.25;
for i=1:params.nWFsToPlot
    plotWaveform(double(squeeze(bckgWF(i,chansToPlot,:)))*params.gain, xPlot, yPlot, 18, 0.2/params.gain, [], 0.75*[1 1 1],p);
    hold on;
end

for i=1:params.nWFsToPlot
    plotWaveform(double(squeeze(theseWF(i,chansToPlot,:)))*params.gain, xPlot, yPlot, 18, 0.2/params.gain, [], neuronColor,p);
    hold on;
end

p.LineWidth = 1; p.alpha = 1;
plotWaveform(medWFuV(:,chansToPlot)', xPlot, yPlot, 18, 0.2/params.gain, [], [0 0 0],p);
box off;
title(sprintf('waveform samples; amp=%.0fµV, SNR=%.2f', wfAmp, snr));

%% plot of PC space
subplot(4,5,[4 5 9 10]);
% figure;
hold off

% h = plot(otherSpikesPCs(otherPCsToPlot,1), otherSpikesPCs(otherPCsToPlot,2), ...
%     '.', 'MarkerSize', 0.05, 'Color', otherColor);
h = plot(otherPCsToPlot(:,1), otherPCsToPlot(:,2), ...
    '.', 'MarkerSize', 3, 'Color', otherColor);
drawnow;
% hMarkers = h.MarkerHandle;
hold on;
% h2 = plot(thesePCs(:,topChans(1)), thesePCs(:,topChans(2)), ...
%     '.', 'MarkerSize', 0.05, 'Color', neuronColor);
h = plot(thesePCsToPlot(:,1), thesePCsToPlot(:,2), ...
    '.', 'MarkerSize', 3, 'Color', neuronColor);
drawnow;
% hMarkers2 = h2.MarkerHandle;

title(sprintf('PC features, iso distance = %.2f, est. contamination = %.2f', stats.isoDistance, stats.mahalContamination))
% set(gca, 'YTickLabel', [], 'XTickLabel', []);
% cEdge = hMarkers.EdgeColorData;
% cEdge(4) = uint8(0.25*255);
% cEdge2 = hMarkers2.EdgeColorData;
% cEdge2(4) = uint8(0.25*255);
% cFace = hMarkers.FaceColorData;
% cFace2 = hMarkers2.FaceColorData;
% addlistener(h,'MarkedClean',...
%     @(ObjH, EventData) keepAlpha(ObjH, EventData, cFace, cEdge));
% addlistener(h2,'MarkedClean',...
%     @(ObjH, EventData) keepAlpha(ObjH, EventData, cFace2, cEdge2));

%% plot of ACG

nSpikes = length(theseST);

subplot(4,5,[14])
hold off;
binSize = 0.02;
acgBins = 1/Fs/2:binSize:0.5;
[n,x] = histdiff(theseST, theseST, acgBins);
n = n./binSize./nSpikes;
stairs([-x(end:-1:1) x]-binSize/2,[n(end:-1:1) n], 'Color', neuronColor, 'LineWidth', 2.0);
yl = ylim(); 
ylim([0 yl(2)]);
title('ACG')
ylabel('spike rate (sp/sec)')
xlabel('time rel. to spike (sec)');

subplot(4,5,[15])
hold off;
binSize = 0.00025;
acgBins = 1/Fs/2:binSize:0.02;
[n,x] = histdiff(theseST, theseST, acgBins);
n = n./binSize./nSpikes;
stairs([-x(end:-1:1) x]-binSize/2,[n(end:-1:1) n], 'Color', neuronColor, 'LineWidth', 2.0);
hold on;
yl2 = ylim();
ylim([0 max([yl(2) yl2(2)])]);
yl2 = ylim();
plot([-0.0015 -0.0015], yl2, 'k--');
plot([0.0015 0.0015], yl2, 'k--');
title(sprintf('ACG zoom, ISI cont = %.2f', stats.isiContamination))
xlabel('time rel. to spike (sec)');

subplot(4,5,[14])
ylim([0 max([yl(2) yl2(2)])]);

%% plot of stability over time
subplot(4,5, 17:19);
hold off;

if isfield(params, 'highlightRegions') && ~isempty(params.highlightRegions)
    % field is nx2    
    for q = 1:size(params.highlightRegions,1)
        hr = params.highlightRegions(1,:);
        fill([hr(1) hr(2) hr(2) hr(1)], [0 0 1000 1000], [0.5 0.5 0.5], 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
        hold on;
    end
end


binSize = 20;
timeBins = 0:binSize:ceil(spikeTimes(end));
[n,x] = hist(theseST, timeBins);
n = n./binSize;

tpc = thesePCs(:,topChans(1));
tpc = (tpc-mean(tpc))./std(tpc)*std(n);
tpc = tpc-min(tpc);

% yyaxis right

h1 = plot(theseST, tpc, '.', 'Color', neuronColor);
hold on;
% ylabel('top PC')

% yyaxis left


h2 = stairs(x,n, 'LineWidth', 2.0, 'Color', 'k');

ymax = max([max(n) max(tpc)]);
ylim([0 ymax]);

title('firing rate over time')
xlabel('time(sec)')
ylabel('firing rate (sp/sec)')

legend([h1 h2], {'top PC', 'firing rate'})

box off;

%% histogram of fit amplitudes 

subplot(4,5,20)

h = histogram(spikeAmps(clu==clusterID), 'Normalization', 'pdf', ...
    'DisplayStyle', 'stairs');
n = get(h, 'Values'); x = get(h, 'BinEdges');
stairs(x(1:end-1), n, 'LineWidth', 2.0);

xlabel('relative spike amplitude')
ylabel('pdf')
title('spike amplitude distribution');


end

% function keepAlpha(src,eventData, FaceColor, EdgeColor)  
%     hm = src.MarkerHandle;
%     hm.EdgeColorData = EdgeColor;
%     hm.FaceColorData = FaceColor;   
% end

