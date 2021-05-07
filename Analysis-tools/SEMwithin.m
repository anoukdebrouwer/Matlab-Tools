function [M,SEM,CI95] = SEMwithin(X)
% SEMwithin Within-subject standard error of the mean
%
% [M,SEM,CI95] = SEMwithin(X) returns the mean M, within-subject standard 
%   error SEM (as described by Cousineau, 2005), and 95% confidence 
%   interval CI95 of the data in X. X should be a matrix of values with a 
%   row for each subject and a column for each condition.
%
% Cousineau, D (2005). Confidence intervals in within-subject designs: A 
%   simpler solution to Loftus and Masson?s method. Tutorial in 
%   Quantitative Methods for Psychology, 1(1):4?45. 
% Url: http://doi.org/10.20982/tqmp.01.1.p042
% 
% See also Sebastiaan Mathot's tutorial here: 
%   http://www.cogsci.nl/blog/tutorials/156-an-easy-way-to-create-graphs-with-within-subject-error-bars

% MIT License
% Copyright (c) 2021 Anouk de Brouwer

% number of subjects
nSubj = size(X,1);

% mean across subjects
M = nanmean(X);

%% Calculate within-subject standard error of mean

% subject mean
subjMean = nanmean(X,2);

% grand mean
grandMean = mean(subjMean);

% correct X
X_corrected = X - repmat(subjMean,1,size(X,2)) + grandMean;

% calculate SEM
SEM = std(X_corrected) / sqrt(nSubj);

% calculate 95% CI
CI95 = SEM * 1.96;
