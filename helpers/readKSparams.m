function params = readKSparams(filename)

fid = fopen(filename);
q = textscan(fid, '%s = %s');
for n = 1:length(q{1})
    try
        params.(q{1}{n}) = eval(q{2}{n});
    catch
        params.(q{1}{n}) = q{2}{n};
    end            
end
