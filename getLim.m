function lim = getLim(X,addPercentage)
% getLim  Range of values for figure axis limits
%
% lim = getLim(X,addPercentage) returns the minimum and maximum value in X,
% minus and plus a percentage (default=5) of the difference between these
% values, to easily set nice figure axis limits.

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

if nargin==1
    addPercentage = 5;
end

minX = nanmin(X(:));
maxX = nanmax(X(:));
rangeX = range(X(:));

lower = minX - rangeX*addPercentage/100;
upper = maxX + rangeX*addPercentage/100;

lim = [lower upper];
