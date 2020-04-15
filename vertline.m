function vertline(x,varargin)
% vertline
% plot vertical line in current figure
%
% vertline(x) plots vertical lines in the current figure at the location 
% specified by x.

% vertline(x,'linespec') sets the line color and style. Default is 'k:'.

% vertline(x,Name,Value) specifies line properties using one or more Name, 
% Value pair arguments. 
% 'color' sets the line color, specified by an RGB triplet or color string.
% 'linestyle' sets the line style. '-' | '--' | ':' | '-.'
% 'linewidth' sets the line width. Default 0.5.

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

% remove NaNs
x = x(~isnan(x));

% parse inputs
p = inputParser;
p.CaseSensitive = false;
addRequired(p,'x',@isnumeric)
addOptional(p,'linespec','',@(x) ischar(x) && length(x)<=3)
addParameter(p,'color','k'); % string or [r g b]
addParameter(p,'linestyle',':',@ischar)
addParameter(p,'linewidth',0.5,@isnumeric)
parse(p,x,varargin{:})

% get values
col = []; 
lin = [];
% if linespec is defined, use value
if ~ismember('linespec',p.UsingDefaults) 
    col = p.Results.linespec(isletter(p.Results.linespec));
    lin = p.Results.linespec(~isletter(p.Results.linespec));
end
% define color if undefined
if isempty(col)
    col = p.Results.color;
end
% define linespec if undefined
if isempty(lin)
    lin = p.Results.linestyle;
end
% define linewidth
width = p.Results.linewidth;

% plot horizontal line(s)
yRange = get(gca,'Ylim');
nLines = length(x);
hold on
for l = 1 : nLines
    plot([x(l) x(l)],yRange,'color',col,'linestyle',lin,'linewidth',width)
end
hold off
