
function [clusterIDs, unitQuality, contaminationRate] = maskedClusterQualityKilosort(resultsDirectory)

try
    pc_features = readNPY([resultsDirectory 'pc_features.npy']); % features of each spike
catch me
    if ~exist([resultsDirectory 'pc_features.npy'])
        fprintf(1, 'PC Features loading failed. File does not exist.\n');
    else
        fprintf(1, 'PC Features loading failed. You may need to clone the npy-matlab repo and add to path.\n');
    end
    rethrow me
end
pc_feature_ind = readNPY([resultsDirectory 'pc_feature_ind.npy']); % the (small) subset of channels corresponding to each template
pc_feature_ind = pc_feature_ind + 1;   % because in Python indexes start at 0

if exist([resultsDirectory 'spike_clusters.npy'])
    spike_clusters = readNPY([resultsDirectory 'spike_clusters.npy']);
    spike_clusters = spike_clusters + 1; % because in Python indexes start at 0
        
    % now we have to construct a new pc_features that has the features
    % arranged according to *cluster* rather than template. 
    spike_templates = readNPY([resultsDirectory 'spike_templates.npy']);
    spike_templates = spike_templates+1;
    
    clusterIDs = unique(spike_clusters);
    nClusters = length(clusterIDs);
    nSpikes = length(spike_clusters);
    nFet = 4; nFetPerChan = size(pc_features,2);
    nTemplates = size(pc_feature_ind,1);
    
    newFet = zeros(nSpikes, nFetPerChan, nFet);
    newFetInds = zeros(nClusters, nFet);
    tempNums = 1:nTemplates;
    
    for c = 1:length(clusterIDs)
        thisID = clusterIDs(c);
        
        theseSpikes = spike_clusters==thisID;
        theseTemplates = spike_templates(theseSpikes);
        inclTemps = unique(theseTemplates); 
        inst = countUnique(theseTemplates); 
        
        thisTemplate = inclTemps(inst==max(inst),1);
        
        theseChans = pc_feature_ind(thisTemplate,1:nFet);
        newFetInds(c,:) = theseChans;
        
        subPCFetInd = pc_features(theseSpikes,:,:);
        
                
        for f = 1:nFet
            thisChanInds = pc_feature_ind==theseChans(f);
            [chanInds,tempsWithThisChan] = find(thisChanInds');
            spikesWithThisFet = ismember(theseTemplates, tempsWithThisChan);
                        
            inclTempsWithThisFet = find(ismember(inclTemps, tempsWithThisChan));
            for t = 1:numel(inclTempsWithThisFet)
                thisSubTemp = inclTemps(inclTempsWithThisFet(t))
                newFet(spikesWithThisFet,:,f) = pc_features(
            
            
        end
    end
else
    fprintf(1, 'warning, spike_clusters does not exist, using spike_templates instead\n');
    spike_clusters = readNPY([resultsDirectory 'spike_templates.npy']); % template # corresponding to each spike
    spike_clusters = spike_clusters + 1; % because in Python indexes start at 0
end



assert(numel(size(pc_features)) == 3)

[clusterIDs, unitQuality, contaminationRate] = maskedClusterQualitySparse(spike_clusters, pc_features, pc_feature_ind);



