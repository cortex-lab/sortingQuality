
function [clusterIDs, unitQuality, contaminationRate] = maskedClusterQualityKilosort(resultsDirectory)

fprintf(1, 'loading data...\n');
%% Precompute the locationsn of files to be loaded
pcFeaturesPath = fullfile(resultsDirectory,'pc_features.npy');
pcFeaturesIndPath = fullfile(resultsDirectory,'pc_feature_ind.npy');
spikeClustersPath = fullfile(resultsDirectory,'spike_clusters.npy');
spikeTemplatesPath = fullfile(resultsDirectory,'spike_templates.npy');


%% Main code.
try
    pc_features = readNPY(pcFeaturesPath); % features of each spike
catch me
    if ~exist(pcFeaturesPath, 'file')
        fprintf(1, 'PC Features loading failed. File does not exist.\n');
    else
        fprintf(1, 'PC Features loading failed. You may need to clone the npy-matlab repo and add to path.\n');
    end
    rethrow me
end
pc_feature_ind = readNPY(pcFeaturesIndPath); % the (small) subset of channels corresponding to each template
pc_feature_ind = pc_feature_ind + 1;   % because in Python indexes start at 0

if exist(spikeClustersPath,'file')
    fprintf('building features matrix from clusters/templates\n')
    spike_clusters = readNPY(spikeClustersPath);
    spike_clusters = spike_clusters + 1; % because in Python indexes start at 0
        
    % now we have to construct a new pc_features that has the features
    % arranged according to *cluster* rather than template. 
    spike_templates = readNPY(spikeTemplatesPath);
    spike_templates = spike_templates+1;
    
    clusterIDs = unique(spike_clusters);
    nClusters = length(clusterIDs);
    nSpikes = length(spike_clusters);
    nFet = 4; nFetPerChan = size(pc_features,2);
    nTemplates = size(pc_feature_ind,1);
    
    newFet = zeros(nSpikes, nFetPerChan, nFet);
    newFetInds = zeros(nClusters, nFet);
    %tempNums = 1:nTemplates;
    
    for c = 1:length(clusterIDs)
%         fprintf(1, '%d/%d\n', c, length(clusterIDs))
        thisID = clusterIDs(c);
        
        theseSpikes = spike_clusters==thisID;
        theseTemplates = spike_templates(theseSpikes);
        [inclTemps, inst] = countUnique(theseTemplates); 
        
        thisTemplate = inclTemps(inst==max(inst),1);
        
        theseChans = pc_feature_ind(thisTemplate,1:nFet);
        
        
        
        newFetInds(c,:) = theseChans;
        
        %subPCFetInd = pc_features(theseSpikes,:,:);
        
                
        for f = 1:nFet
            thisChanInds = pc_feature_ind==theseChans(f);
            [chanInds,tempsWithThisChan] = find(thisChanInds');
            %spikesWithThisFet = ismember(theseTemplates, tempsWithThisChan);
                        
            inclTempsWithThisFet = find(ismember(inclTemps, tempsWithThisChan));
            for t = 1:numel(inclTempsWithThisFet)
                thisSubTemp = inclTemps(inclTempsWithThisFet(t));
                thisTfetInd = chanInds(tempsWithThisChan==thisSubTemp);
                newFet(theseSpikes&spike_templates==thisSubTemp,:,f) = ...
                    pc_features(theseSpikes&spike_templates==thisSubTemp,:,thisTfetInd);
            end
            
            
        end
    end
    
    pc_features = newFet;
    pc_feature_ind = newFetInds;
else
    fprintf(1, 'warning, spike_clusters does not exist, using spike_templates instead\n');
    spike_clusters = readNPY(spikeTemplatesPath); % template # corresponding to each spike
    spike_clusters = spike_clusters + 1; % because in Python indexes start at 0
end



assert(numel(size(pc_features)) == 3)

fprintf(1, 'computing cluster qualities...\n');
[clusterIDs, unitQuality, contaminationRate] = maskedClusterQualitySparse(spike_clusters, pc_features, pc_feature_ind);



