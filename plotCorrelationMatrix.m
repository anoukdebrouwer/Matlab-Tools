function [R,P] = plotCorrelationMatrix(X,labels)
% plotCorrelationMatrix
% Plot correlation matrix
%
% plotCorrelationMatrix(X,labels) creates a plot of the matrix of Pearson
% correlation coefficients between the variables in X.
% plotCorrelationMatrix displays the correlations between each two
% variables as a colored cell in a table, ranging from red for positive
% correlations to white for zero correlations, to blue for negative
% correlations. The correlation coefficient is displayed in each cell with
% asterisks indicating the p value. Only the below-diagonal correlations
% are displayed for simplicity.
%
% [R,P] = plotCorrelationMatrix(X,labels) also returns the matrix of
% correlation coefficients and the matrix of p-values for testing the
% hypothesis that there is no relationship between the variables (null
% hypothesis).

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

cla reset % clear and reset figure axes

nVar = size(X,2);

% if labels are not provided, use numbers
if nargin==1
    labels = num2cell(1:nVar);
end

% compute correlations
[R,P] = corrcoef(X,'rows','pairwise');

% hide above-diagonal values by making r values above the diagonal 0 and
% p values above the diagonal 1
rowIndex = repmat(1:nVar,1,nVar);
colIndex = repelem(1:nVar,nVar);
R(colIndex>rowIndex) = 0;
P(colIndex>rowIndex) = 1;

% create colormap
blueWhiteRed = [createColorGradient([0 0 1],[1 1 1],32);...
    createColorGradient([1 1 1],[1 0 0],32)];

% create image
imagesc(R)
colormap(blueWhiteRed)
caxis([-1,1])
c = colorbar('eastoutside');
c.Label.String = 'r value';
set(gcf,'color','w')
set(gca,'TickLabelInterpreter','none')
set(gca,'XTick',1:nVar)
set(gca,'XTickLabels',labels)
nChar = cellfun(@length,labels);
if any(nChar>1)
    set(gca,'XTickLabelRotation',45)
end
set(gca,'YTick',1:nVar)
set(gca,'YTicklabels',labels)
set(gca,'TickLength',[0 0])
axis square; box off
title('Correlations')

% get estimate of text size displaying r-values and asterisks for
% significance and determine the maximum number of variables for which
% text will fit
t = text(1,1,' 0.00*** ');
t.Units = 'normalized';
textSize = t.Extent;
delete(t);
maxVar(1) = floor(1/textSize(3));
t = text(1,1,' *** ');
t.Units = 'normalized';
textSize = t.Extent;
maxVar(2) = floor(1/textSize(3));
delete(t);
% add text if it fits
if nVar<=maxVar(2)
    for row = 1 : nVar
        for col = 1 : nVar
            if row>col % below diagonal
                if nVar<=maxVar(1)
                    str = num2str(R(row,col),'%.2f');
                else
                    str = '';
                end;
                if P(row,col)<0.001
                    str = [str '***'];
                elseif P(row,col)<0.01
                    str = [str '**'];
                elseif P(row,col)<0.05
                    str = [str '*'];
                end
                text(col,row,str,'HorizontalAlignment','center');
            end
        end
    end
end
% add legend for significance
if nVar<=maxVar(1)
    text(nVar*0.8,1,{'*   p<0.05','**  p<0.01','*** p<0.001'},...
        'HorizontalAlignment','Left')
end

end

function colorGradient = createColorGradient(startColor,endColor,nColors)

colorGradient = zeros(nColors,3);
for i = 1 : 3
    if endColor(i)~=startColor(i)
        colorGradient(:,i) = startColor(i) : (endColor(i)-startColor(i))/(nColors-1) : endColor(i);
    elseif endColor(i)==startColor(i)
        colorGradient(:,i) = startColor(i);
    end
end

end