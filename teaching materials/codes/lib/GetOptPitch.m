function pitch = GetOptPitch(WSPfilt)

persistent wpdata
if isempty(wpdata)
    formatSpec = '%f';
    wpdata = dlmread('wpdata.100', formatSpec, 1, 0);
end

pitch = interp1(wpdata(:,1),wpdata(:,2),WSPfilt,'linear');