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
        [tr,day] = find(outlier);
        % plot reported values and removed outliers
        figure(fig1); clf
        plot([Results.explicitAngle_outliersRemoved],'o'); hold on
        for i = 1 : length(tr)
            plot(tr(i),Results.explicitAngle(tr(i),day(i)),'kx');
        end
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
    bin1 = pearlNumAngle(2,2);
    if vmr>0; bin1 = -bin1; end
    firstAimFixBin(s) = nanmax([find(aimFixAngle(firstRotBin:end,s,1)<=bin1,1) NaN]);
    firstAimReportBin(s) = nanmax([find(reportAngle(firstRotBin:end,s,1)<=bin1,1) NaN]);
    
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
        'cursorRotation','meanOrMedian','RT','handAngle','aimFixAngle','aimFixAngle_opp',...
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

% colors
colors = brewermap(10,'Paired');
colors = colors([1,2,5,6,3,4,9,10],:); % blue, red, green, purple

% plot
scaledFigure(nDays*0.75,0.5*4);
subplotNo = 1;
for i = 1 : 4
    for d = 1 : nDays
        subplot(4,nDays,subplotNo); hold on
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
        if i==1
            plot(iMidBin,handAngle(:,:,d),'.-','color',colors(1,:));
            plot(iMidBin,nanmean(handAngle(:,:,d),2),'.-','color',colors(2,:),'linewidth',2);
            figTitle = 'Hand angle';
        elseif i==2
            plot(iMidBin,reportAngle(:,:,d),'.-','color',colors(3,:));
            plot(iMidBin,nanmean(reportAngle(:,:,d),2),'.-','color',colors(4,:),'linewidth',2);
            figTitle = 'Report angle';
        elseif i==3
            plot(iMidBin,aimFixAngle(:,:,d),'.-','color',colors(7,:));
            plot(iMidBin,aimFixAngle_opp(:,:,d),'.-','color',colors(7,:));
            mAimFixAngle = nanmean([aimFixAngle(:,:,d) aimFixAngle_opp(:,:,d)],2);
            plot(iMidBin,mAimFixAngle,'.-','color',colors(8,:),'linewidth',2);
            nAimFixators = sum(~isnan(aimFixAngle(:,:,d)),2);
            str = ['max n = ' num2str(max(nAimFixators))];
            text(iMidBin(firstRotBin),0.5*vmr,str)
            figTitle = 'Aimpoint fixation angle';
        elseif i==4
            % from report
            plot(iMidBin,implicitAngle(:,:,d),'.-','color',colors(5,:));
            mImplAngle = nanmean(implicitAngle(:,:,d),2);
            pr = plot(iMidBin,mImplAngle,'.-','color',colors(6,:),'linewidth',2);
            % from fixations
            implicitAngle_fromFix(~isnan(mImplAngle),:,d) = NaN;
            plot(iMidBin,implicitAngle_fromFix(:,:,d),'.--','color',colors(5,:));
            mImplAngle_fromFix = nanmean(implicitAngle_fromFix(:,:,d),2);
            pf = plot(iMidBin,mImplAngle_fromFix,'.--','color',colors(6,:),'linewidth',2);
            figTitle = 'Implicit angle';
            if d==1 && sum(~isnan(mImplAngle))>1 && sum(~isnan(mImplAngle_fromFix))>1
                legend([pr pf],{'from report','from fixation'},'Location','best','Orientation','horizontal');
            end
        end
        % axes and lables
        hold off
        vertline(breakTrials,':','color',[0.7 0.7 0.7])
        vertline(iNewBlock,'k:')
        if i==4
            xlabel('Trial')
        end
        if d==1
            ylabel('Angle (deg)')
        end
        if nDays==1
            title(figTitle)
        elseif nDays==2
            %title(['Day ' num2str(d) ' - ' figTitle])
            title([figTitle ' [day ' num2str(d) ']' ])
        end
        subplotNo = subplotNo+1;
    end
end
figTitle = [expFolder.name ' - group results (n=' num2str(nSubj) ')'];
suplabel(figTitle,'t');

% save figure
figName = [expFolder.name '_binnedAngles'];
if exist([saveToPath figName],'file') == 2 % check if file does not exist yet
    disp(['A figure named ' figName ' already exists.'])
    overwrite = input('Do you want to overwrite it? Yes(1) or no(0): ');
else
    overwrite = 1;
end
if overwrite == 1
    saveFigAsPDF(figName,12)
end
