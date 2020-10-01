function saveCopyOfCode(mFilePath,saveToPath)
% saveCopyOfCode Save copy of Matlab code used to process data

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

% get program dependencies
[fList,~] = matlab.codetools.requiredFilesAndProducts([mFilePath '.m']);
fList = [[mFilePath '.m'] fList];
fList = fList(cellfun(@isempty,strfind(fList,'saveCopyOfCode'))); % discard this function
fList = fList(cellfun(@isempty,strfind(fList,'/Plotting-tools/'))); % discard plotting tools
fList = fList(cellfun(@isempty,strfind(fList,'/Analysis-tools/'))); % discard analysis tools
fList = unique(fList); % remove duplicates

% loop over files and save a copy with the today's date in the name
fileDate = datestr(now,'yyyymmmmdd');
for f = 1 : length(fList)
    mFilePath = fList{f}; % overwrite file path with current file
    iSlash = strfind(mFilePath,'/');
    mfileName = mFilePath(iSlash(end)+1:end);
    mfileCopyName = sprintf('%s_copy_%s',mfileName,fileDate);
    copyfile(mFilePath,[saveToPath mfileCopyName]);
end