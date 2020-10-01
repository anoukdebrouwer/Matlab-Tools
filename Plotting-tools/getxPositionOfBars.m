function [xBars] = getxPositionOfBars(barHandle)
% getxPositionOfBars
% Get the x position of bars in a bar plot.

% [xBars] = getxPositionOfBars(barHandle) returns an array with the x 
% positions of the centers of the bars in a bar plot, e.g. for plotting 
% error bars.

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

xData = get(barHandle,'XData');
if iscell(xData)
    xData = cell2mat(xData)';
end
xOffset = get(barHandle,'XOffset');
if iscell(xOffset)
    xOffset = cell2mat(xOffset)';
end
xBars = xData + repmat(xOffset,size(xData,1),1);
