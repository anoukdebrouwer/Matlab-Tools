function folderName = getFolderName(path,level)
% getFolderName Get name of folder from path.
%
% getFolderName(path) returns the name of the bottom folder in path.
%
% getFolderName(path,level) gets the name of the folder at a
% specified level above the bottom folder.

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

if nargin==1
    level = 0;
end

path_split = strsplit(path,filesep); % split path
path_split = path_split(cellfun(@(x) ~isempty(x),path_split)); % remove empty strings
folderName = path_split{end-level}; % get folder name