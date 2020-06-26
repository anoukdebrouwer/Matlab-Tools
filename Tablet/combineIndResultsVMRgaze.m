function combineIndResultsVMRgaze(meanOrMedian)
% combineIndResultsVMRgaze Combine individual binned data of VMR gaze
% experiments into a group data file.
%
% combineResultsVMRgaze loads the individual data files and creates and
% saves matrices with the binned hand angle, reported aiming angle, aim
% fixation angle, and implicit angle for all participants. It also creates
% and saves a figure with the individual learning curves (hand angles).
%
% For 2-day experiments, it is assumed that the experimental blocks on both
% days are identical (although day 2 doesn't need to include all day-1 blocks).

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

if nargin==0
    meanOrMedian = 1;
end

% select experiment data folder
projectPath = '/Users/anouk/Documents/ExpsTablet/';
expFolder = selectFiles(projectPath,'folders');
while ~exist([projectPath expFolder.name '/1_RawData/'],'dir');
    projectPath = [projectPath expFolder.name '/'];
    expFolder = selectFiles(projectPath,'folders'); % look in subfolders
end
dataPath = [projectPath expFolder.name '/3_Results/'];
saveToPath = [projectPath expFolder.name '/3_Results/'];
cd(dataPath)

% select subject data files
subjFiles = selectFiles('*.mat');
nSubj = length(subjFiles);
subj = cell(1,nSubj);

% load experiment details
detailsFile = dir([projectPath expFolder.name '/ExpDetails*.mat']);
load([projectPath expFolder.name '/' detailsFile.name])
nBins = nTrials/nTargets;
vmr = visuomotorRotation;

% open figures
fig1 = figure;

% preallocate
RT                      = NaN(nBins,nSubj,nDays);
handAngle               = NaN(nBins,nSubj,nDays);
reportAngle             = NaN(nBins,nSubj,nDays);
implicitAngle           = NaN(nBins,nSubj,nDays);
implicitAngle_fromFix   = NaN(nBins,nSubj,nDays);
aimFixAngle             = NaN(nBins,nSubj,nDays);
aimFixAngle_opp         = NaN(nBins,nSubj,nDays);
firstAimReportBin       = NaN(nSubj,1);
firstAimFixBin          = NaN(nSubj,1);

%% Loop over subjects

for s = 1 : nSubj
    
    load(subjFiles(s).name)
    disp(['Loaded ' subjFiles(s).name])
    subj{s} = ExpDetails(1).subjName;
    
    % get general experiment information
    if s==1
        iBin = Results.trialNo_bins(:,end);
        blockNo = Results.blockNo_bins(:,1);
        cursorRotation = Results.cursorRotation(iBin);
        rotBins = cursorRotation ~= 0;
        woBlock = find(ismember(lower(blockNames),'washout'));
        woBins = blockNo==woBlock;
    end
    
    % check if means and medians of bins were calculated
    % (older data only has means)
    if ~isfield(Results,'mnBinnedRT')
        meanOrMedian = 1;
        Results.mnBinnedRT = Results.binnedRT;
        Results.mnBinnedHitAngle_hand = Results.binnedHitAngle_hand;
        Results.mnBinnedFixAngle_closestA = Results.binnedFixAngle_closestA;
        Results.mnBinnedFixAngle_closestOppA = Results.binnedFixAngle_closestOppA;
        Results.mnBinnedExplicitAngle = Results.binnedExplicitAngle;
        Results.mnBinnedExplicitAngle_outliersRemoved = Results.binnedExplicitAngle_outliersRemoved;
    end
    
    %% Fill matrix with binned data
    
    if meanOrMedian==1 % get means
        RT(:,s,:) = Results.mnBinnedRT;
        handAngle(:,s,:) = Results.mnBinnedHitAngle_hand;
        if any(Results.analyzedGazeData(:))
            aimFixAngle(rotBins,s,:) = Results.mnBinnedFixAngle_closestA(rotBins,:);
            aimFixAngle_opp(woBins,s,:) = Results.mnBinnedFixAngle_closestOppA(woBins,:);
        end
    else % get medians
        RT(:,s,:) = Results.medBinnedRT;
        handAngle(:,s,:) = Results.medBinnedHitAngle_hand;
        if any(Results.analyzedGazeData(:))
            aimFixAngle(rotBins,s,:) = Results.medBinnedFixAngle_closestA(rotBins,:);
            aimFixAngle_opp(woBins,s,:) = Results.medBinnedFixAngle_closestOppA(woBins,:);
        end
    end
    
    % check removed outliers in reported aiming angles
    keepOutliers = false;
    if meanOrMedian==1
        outlierDiff = ~isnan(Results.mnBinnedExplicitAngle) & ...
            (Results.mnBinnedExplicitAngle ~= Results.mnBinnedExplicitAngle_outliersRemoved);
    else
        outlierDiff = ~isnan(Results.medBinnedexplicitAngle) & ...
            (Results.medBinnedExplicitAngle ~= Results.medBinnedExplicitAngle_outliersRemoved);
    end
    if any(outlierDiff(:));
        outlier = ~isnan(Results.explicitAngle) & ...
            (Results.explicitAngle ~= Results.explicitAngle_outliersRemoved);
        [i,j] = find(outlier);
        % plot reported values and removed outliers
        figure(fig1); clf
        plot([Results.explicitAngle_outliersRemoved],'o'); hold on
        plot(i,Results.explicitAngle(i,j),'kx');
        horline([-vmr,0,vmr])
        if nDays==2
            legend('Day1','Day2','outlier')
        end
        set(gca,'ytick',-180:45:180)
        xlabel('Trial'); ylabel('Angle (deg)')
        title('Reported aiming angles and removed outliers');
        keepOutliers = input('Discard(0) or keep(1) outliers? ');
    end
    if keepOutliers
        if meanOrMedian==1
            reportAngle(:,s,:) = Results.mnBinnedExplicitAngle;
        else
            reportAngle(:,s,:) = Results.medBinnedExplicitAngle;
        end
        keyboard
    else
        if meanOrMedian==1
            reportAngle(:,s,:) = Results.mnBinnedExplicitAngle_outliersRemoved;
        else
            reportAngle(:,s,:) = Results.medBinnedExplicitAngle_outliersRemoved;
        end
    end
    
    % calculate implicit angle
    implicitAngle(:,s,:) = handAngle(:,s,:)-reportAngle(:,s,:);
    implicitAngle_fromFix(rotBins,s,:) = handAngle(rotBins,s,:)-aimFixAngle(rotBins,s,:);
    implicitAngle_fromFix(woBins,s,:) = handAngle(woBins,s,:)-aimFixAngle_opp(woBins,s,:);
    
    % find first bin in which an aiming strategy was present
    firstRotBin = find(rotBins,1);
    firstAimFixBin(s) = nanmax([find(aimFixAngle(firstRotBin:end,s,1)<0,1) NaN]);
    firstAimReportBin(s) = nanmax([find(reportAngle(firstRotBin:end,s,1)<0,1) NaN]);
    
end % end of loop over subjects

%% Save data file

fileName = ['Results_binnedAngles_' expFolder.name '.mat'];

% check if file does not exist yet
if exist([saveToPath fileName],'file') == 2
    disp(['A file named ' fileName ' already exists.'])
    overwrite = input('Do you want to overwrite it? Yes(1) or no(0): ');
else
    overwrite = 1;
end

% save file
if overwrite == 1
    save([saveToPath fileName],'subj','blockNames','blockNo','iBin','nTrials',...
        'cursorRotation','meanOrMedian','RT','handAngle','aimFixAngle',...
        'reportAngle','implicitAngle','implicitAngle_fromFix','firstAimReportBin','firstAimFixBin');
    disp(['Saved ' saveToPath fileName])
    
    % save copy of Matlab code used to process data
    mFilePath = mfilename('fullpath');
    saveCopyOfCode(mFilePath,saveToPath)
else
    disp('Results not been saved')
end

%% Plot individual learning curves

iMidBin = iBin+3;
iRotationOnOff = find(diff(Results.cursorRotation(:,1))) + 0.5;
if mod(length(iRotationOnOff),2)
    iRotationOnOff = [iRotationOnOff; nTrials+0.5];
end
iNewBlock = find(diff(Results.blockNo(:,1))) + 0.5;

% plot
scaledFigure(1.5,0.75*nDays);
for d = 1 : nDays
    subplot(nDays,1,d); hold on
    xlim([0 iMidBin(end)+1]); set(gca,'XTick',0:80:nTrials)
    ylim([-2*vmr +vmr]); yl = ylim; set(gca,'Ytick',-180:45:180)
    % color background when rotation is on, works in R2015b but not in 2017a
    for r = 1 : length(iRotationOnOff)/2
        x = iRotationOnOff(r*2-1:r*2);
        a = area([x; flipud(x)],[yl(1) yl(1) yl(2) yl(2)],'LineStyle','none');
        a.FaceColor = [0.95 0.95 0.95];
    end
    plot([0 nTrials],[0 0],'k'); % draw line at zero
    plot([0 nTrials],[-vmr -vmr],'k'); % draw line at hand target angle
    % plot learning curves
    pt = plot(iMidBin,handAngle(:,:,d));
    % axes and lables
    hold off
    vertline(breakTrials,'g:')
    vertline(iNewBlock,'k:')
    xlabel('Trial')
    ylabel('Angle (deg)')
    if nDays==2
        if d==1
            title('Day 1')
        elseif d==2
            title('Day 2')
        end
    end
end
figTitle = [expFolder.name ' - Hand angles of all partipants (n=' num2str(nSubj) ')'];
if nDays==2
    suplabel(figTitle,'t');
else
    title(figTitle,'interpreter','none');
end

% save figure
figName = [expFolder.name '_binnedHandAngles'];
if exist([saveToPath figName],'file') == 2 % check if file does not exist yet
    disp(['A figure named ' figName ' already exists.'])
    overwrite = input('Do you want to overwrite it? Yes(1) or no(0): ');
else
    overwrite = 1;
end
if overwrite == 1
    saveFigAsPDF(figName,12)
end
