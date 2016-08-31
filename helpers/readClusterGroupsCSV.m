

function [cids, cgs] = readClusterGroupsCSV(filename)
%function [cids, cgs] = readClusterGroupsCSV(filename)

fid = fopen(filename);
C = textscan(fid, '%s%s');
fclose(fid);

cids = cellfun(@str2num, C{1}(2:end), 'uni', false);
ise = cellfun(@isempty, cids);
cids = [cids{~ise}];

isUns = cellfun(@(x)strcmp(x,'unsorted'),C{2}(2:end));
isMUA = cellfun(@(x)strcmp(x,'mua'),C{2}(2:end));
isGood = cellfun(@(x)strcmp(x,'good'),C{2}(2:end));
cgs = zeros(size(cids));

cgs(isMUA) = 1;
cgs(isGood) = 2;
cgs(isUns) = 3;