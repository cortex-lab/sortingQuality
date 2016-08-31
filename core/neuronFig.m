

function f = neuronFig(clusterID, spikeTimes, clu, sparsePCfeat, stats, params)
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
extractST = round(theseST(randperm(length(theseST), params.nWFsToLoad))*params.Fs); 
extractBckg = randperm(nSamp-numel(wfWin)-2, params.nWFsToPlot)-wfWin(1); % *samples* of background spikes

mmf = memmapfile(params.filename, 'Format', {params.dataType, [nChInFile nSamp], 'x'});

theseWF = zeros(params.nWFsToLoad, nCh, nWFsamps);
bckgWF = zeros(params.nWFsToPlot, nCh, nWFsamps);
for i=1:params.nWFsToLoad
     tempWF = mmf.Data.x(1:nChInFile,extractST(i)+wfWin(1):extractST(i)+wfWin(end));
     theseWF(i,:,:) = tempWF(params.chanMap+1,:);
end

for i = 1:params.nWFsToPlot %
     tempWF = mmf.Data.x(1:nChInFile,extractBckg(i)+wfWin(1):extractBckg(i)+wfWin(end));
     bckgWF(i,:,:) = tempWF(params.chanMap+1,:);
     
end

medWF = squeeze(median(double(theseWF),1))';
medWFuV = medWF.*params.gain;


%%

% choose channels to plot

chanAmps = max(medWFuV)-min(medWFuV);
maxChan = find(chanAmps==max(chanAmps),1);
maxXC = xc(maxChan); maxYC = yc(maxChan);
chanDistances = ((xc-maxXC).^2 + (yc-maxYC).^2).^0.5;
chansToPlot = chanDistances<params.plotDistance;

wfAmp = max(chanAmps);

bckgMaxChan = double(squeeze(bckgWF(:,maxChan,:)));
snr = wfAmp./std(bckgMaxChan(:));

%% compute PC stuff

thesePCs = sparsePCfeat(clu==clusterID,:);
meanPC = mean(thesePCs);
[~, ii] = sort(meanPC, 2, 'descend');

topChans = ii(1:2);
otherSpikesIncl = sparsePCfeat(:,topChans(1))~=0 & sparsePCfeat(:,topChans(2))~=0;
otherSpikesPCs = sparsePCfeat(otherSpikesIncl, topChans);

%%
f = figure;
set(f, 'Color', 'w');

% neuronColor = [0.8500    0.3250    0.0980]; % the orange you get with second element in default color order
neuronColor = [ 0    0.4470    0.7410]; % blue that's first default
otherColor = [0.4660    0.6740    0.1880];

% plot of location on probe
subplot(4,5,[1 6 11 16])
plotAsProbe(-chanAmps, xc, yc, colormap(gray), 16, 40)
hold on;
minx = min(xc)-8; maxx = max(xc)+8; miny = min(yc)-20; maxy = max(yc)+20;
plot([minx maxx maxx minx minx], [maxy maxy miny miny maxy], 'k');
axis off
title('neuron position')
% xlabel('space (µm)')
% ylabel('space (µm)')

% plot of waveforms, zoomed around location
subplot(4,5, [2 3 7 8 12 13]);
xPlot = xc(chansToPlot);
yPlot = yc(chansToPlot);
p.LineWidth = 0.01; p.alpha = 0.25;
for i=1:params.nWFsToPlot
    plotWaveform(double(squeeze(bckgWF(i,chansToPlot,:))), xPlot, yPlot, 18, 0.2, [], 0.75*[1 1 1],p);
    hold on;
end

for i=1:params.nWFsToPlot
    plotWaveform(double(squeeze(theseWF(i,chansToPlot,:))), xPlot, yPlot, 18, 0.2, [], neuronColor,p);
    hold on;
end

p.LineWidth = 1; p.alpha = 1;
plotWaveform(medWF(:,chansToPlot)', xPlot, yPlot, 18, 0.2, [], [0 0 0],p);
box off;
title(sprintf('waveform samples; amp=%.0fµV, SNR=%.2f', wfAmp, snr));

%% plot of PC space
subplot(4,5,[4 5 9 10]);
% figure;
hold off
otherPCsToPlot = randperm(size(otherSpikesPCs,1), params.nPCsToPlot);
h = plot(otherSpikesPCs(otherPCsToPlot,1), otherSpikesPCs(otherPCsToPlot,2), ...
    '.', 'MarkerSize', 0.05, 'Color', otherColor);
drawnow;
% hMarkers = h.MarkerHandle;
hold on;
h2 = plot(thesePCs(:,topChans(1)), thesePCs(:,topChans(2)), ...
    '.', 'MarkerSize', 0.05, 'Color', neuronColor);
drawnow;
% hMarkers2 = h2.MarkerHandle;

title(sprintf('PC features, iso distance = %.2f', stats.isoDistance))

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
plot([-0.002 -0.002], yl2, 'k--');
plot([0.002 0.002], yl2, 'k--');
title(sprintf('ACG zoom, ISI cont = %.2f', stats.isiContamination))
xlabel('time rel. to spike (sec)');

subplot(4,5,[14])
ylim([0 max([yl(2) yl2(2)])]);

%% plot of stability over time
subplot(4,5, 17:20);
hold off;

binSize = 20;
timeBins = 0:binSize:ceil(spikeTimes(end));
[n,x] = hist(theseST, timeBins);
n = n./binSize;

tpc = thesePCs(:,topChans(1));
tpc = (tpc-mean(tpc))./std(tpc)*std(n);
tpc = tpc-min(tpc);

% yyaxis right

plot(theseST, tpc, '.', 'Color', neuronColor);
hold on;
% ylabel('top PC')

% yyaxis left


stairs(x,n, 'LineWidth', 2.0, 'Color', 'k');
title('firing rate over time')
xlabel('time(sec)')
ylabel('firing rate (sp/sec)')

legend({'top PC', 'firing rate'})

box off;

end

% function keepAlpha(src,eventData, FaceColor, EdgeColor)  
%     hm = src.MarkerHandle;
%     hm.EdgeColorData = EdgeColor;
%     hm.FaceColorData = FaceColor;   
% end

