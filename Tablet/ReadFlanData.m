function [Data]=ReadFlanData(FileName)
% ReadFlanData  Text import function for Flanagan lab data

% Jason Gallivan

% open file
fid = fopen(FileName);

% read number of header lines from first row
C=textscan(fid,'%s %d',1);
nHeaderLines = C{2};

% load data into workspace
Data=importdata(FileName,'\t',nHeaderLines+1);

% close file
fclose(fid);