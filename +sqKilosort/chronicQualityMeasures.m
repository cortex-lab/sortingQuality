
function recStats = chronicQualityMeasures(ksDir, rawFilename, inclChans, gain)
% function recStats = chronicQualityMeasures(ksDir, rawFilename, inclChans, gain)
% Computes the agreed-upon quality statistics for showing probe performance
% over time for chronically implanted probes:
% 1) Rate of spikes >50uV
% 2) Rate of spikes >6*MAD
% 3) spike SNR, the peak amplitude divided by Vrms
%
% Inputs:
% - ksDir, a path to the folder with kilosort outputs
% - rawFilename, a path to the raw data file
% - inclChans, which channels should be part of the counting (1:374 will
% include all for a Neuropixels probe)
% - gain, the gain setting on the probe

%%
% ksDir = 'J:\M160731_MOEC\2016-08-16';
% rawFilename = 'Z:\multichanspikes\M160731_MOEC\2016-08-16\M160731_MOEC_2016-08-16_g0_t0.imec.ap_CAR.bin'; 
% inclChans = 1:330;
% gain = 500;

gainFactor = 0.6/512/gain*1e6;

%% step 1: compute MAD for each channel

pars = loadParamsPy(fullfile(ksDir, 'params.py'));
Fs = pars.sample_rate;
nCh = pars.n_channels_dat;

% read 10 sec of data
fid = fopen(rawFilename, 'r');
rawDat = fread(fid, [nCh Fs*10], 'int16=>double');
fclose(fid);

chanMap = readNPY(fullfile(ksDir, 'channel_map.npy'));

% scale by gain
rawDat = rawDat(chanMap+1,:).*gainFactor;

% compute MAD
madPerChannel = mad(rawDat',1)';
recStats.madPerChannel = madPerChannel;
% note that because the lsb of the 10-bit adc is 2.34µV at gain=500, this
% number will come out at only a few different values across all channels
% (i.e. 9.4, 11.7, 14.1, etc). Similar argument will apply for probes other
% than neuropixels. 

%% step 2: compute peak channel and amplitude of every spike

spikeTemplates = readNPY(fullfile(ksDir, 'spike_templates.npy')); % note: zero-indexed
tempScalingAmps = readNPY(fullfile(ksDir, 'amplitudes.npy'));

if exist(fullfile(ksDir, 'spike_clusters.npy'), 'file')
    clu = readNPY(fullfile(ksDir, 'spike_clusters.npy'));
else
    clu = spikeTemplates;
end

if exist(fullfile(ksDir, 'cluster_groups.csv'), 'file') 
    [cids, cgs] = readClusterGroupsCSV(fullfile(ksDir, 'cluster_groups.csv'));
    
    noiseClusters = cids(cgs==0);
    
    spikeTemplates = spikeTemplates(~ismember(clu, noiseClusters));
    tempScalingAmps = tempScalingAmps(~ismember(clu, noiseClusters));
       
end

temps = readNPY(fullfile(ksDir, 'templates.npy'));

winv = readNPY(fullfile(ksDir, 'whitening_mat_inv.npy'));

% unwhiten the templates
tempsUnW = zeros(size(temps));
for t = 1:size(temps,1)
    tempsUnW(t,:,:) = squeeze(temps(t,:,:))*winv;
end

% The amplitude on each channel is the negative peak
tempChanAmps = -squeeze(min(tempsUnW,[],2));

% add epsilon to every channel amp to guarantee there is just one max
% channel for each template (break ties)
tempChanAmps = tempChanAmps+rand(size(tempChanAmps))/1e9; 

% The template amplitude is the amplitude of its largest channel
tempAmps = max(tempChanAmps,[],2);

% The template's assigned channel is the one with the largest amplitude
tempMaxes = bsxfun(@eq, tempChanAmps, tempAmps);
[tempPeakChannel,~] = find(tempMaxes');

% now compute the amplitude for each spike
spikeAmps = tempAmps(spikeTemplates+1).*tempScalingAmps.*gainFactor; 
% and peak channel for each spike
spikePeakChannel = tempPeakChannel(spikeTemplates+1); 

%recStats.spikeAmps = spikeAmps;
%recStats.spikePeakChannel = spikePeakChannel;

%% compile three metrics 1) rate of spikes >50µV

d = dir(rawFilename); nSamp = d.bytes/2/385; recDuration = nSamp/Fs;

[chanNums, countsPerChan] = countUnique(spikePeakChannel(spikeAmps>50)); 
spikeCounts = zeros(size(chanMap)); 
spikeCounts(chanNums) = countsPerChan;
spikeRates50 = spikeCounts./recDuration;

recStats.spikeRates50 = spikeRates50;
recStats.spikeRates50iqr = prctile(spikeRates50(inclChans), [25 50 75]);
recStats.spikeRates50total = sum(spikeRates50(inclChans));
recStats.recDuration = recDuration;

%% compile three metrics 2) rate of spikes > 6*MAD

% get the MAD corresponding to the peak channel of the spike
spikeMAD = madPerChannel(spikePeakChannel); 

spikeSNR = spikeAmps./spikeMAD; 

[chanNums, countsPerChan] = countUnique(spikePeakChannel(spikeSNR>6)); 
spikeCounts = zeros(size(chanMap)); 
spikeCounts(chanNums) = countsPerChan;
spikeRates6MAD = spikeCounts./recDuration;

recStats.spikeRates6MAD = spikeRates6MAD;
recStats.spikeRates6MADiqr = prctile(spikeRates6MAD(inclChans), [25 50 75]);
recStats.spikeRates6MADtotal = sum(spikeRates6MAD(inclChans));

%% compile three metrics 3) SNR of spikes relative to MAD
% From James: To calculate SNR, I define Vp as the peak negative voltage 
% (absolute value) of the largest spike waveform per spiking event. I then 
% divide this by the RMS voltage (Vrms=MAD/0.6745) ... I do not count a 
% given spike more than once and I only use the largest spike per given 
% spiking event. I pool all spiking events detected from all the sites, 
% and compute the median and IQR across the pooled population.

% Also, since james uses a 50uV detection threshold, and this analysis
% depends on how many spikes you detect, we need to exclude spikes with amp
% < 50uV

inclSpikes = spikeAmps>50; 
spikeSNR = spikeAmps(inclSpikes)./(spikeMAD(inclSpikes)/0.6745);

recStats.spikeSNRiqr = prctile(spikeSNR, [25 50 75]);
%recStats.spikeSNR = spikeSNR;

return;


%% test: do amplitudes computed from templates match the raw data?

ss = readNPY(fullfile(ksDir, 'spike_times.npy'));

ss = ss(ss<Fs*10); % spikes for which we already have the raw data
thesePeakChans = spikePeakChannel(ss<Fs*10); 

% so for each spike we want to extract the minimum raw data value in a
% small window around the spike time on its peak channel
rawAmps = zeros(size(ss));
for q = 1:length(ss)
    rawAmps(q) = -min(rawDat(thesePeakChans(q), max(1,ss(q)-30):min(Fs*10,ss(q)+30)));
end

figure; 
plot(rawAmps, spikeAmps(ss<Fs*10), '.');
hold on; plot([0 max(rawAmps)], [0 max(rawAmps)], 'k--');
xlabel('raw empirical amplitude'); 
ylabel('amplitude from template');
% looks good. consistently a very slight underestimate because taking the
% peak of a noisy signal will overestimate the true peak (i.e. the raw data
% is an overestimate of the true peak that's partly corrected by using the
% templates)

