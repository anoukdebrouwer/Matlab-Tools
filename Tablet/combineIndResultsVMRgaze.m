function combineIndResultsVMRgaze(meanOrMedian,createPlots)
% combineIndResultsVMRgaze Bin individual data of VMR gaze experiments and
% combine into a group data file.
%
% combineResultsVMRgaze loads the individual data files and creates and
% saves matrices with the binned hand angle, reported aiming angle, aim
% fixation angle, and implicit angle for all participants.
%
% For 2-day experiments, it is assumed that the experimental blocks on both
% days are identical (although day 2 doesn't need to include all day-1 blocks).

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

if nargin==0
    meanOrMedian = 1;
    createPlots = false;
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
sdRT                    = NaN(nBins,nSubj,nDays);
cursorAngle             = NaN(nBins,nSubj,nDays);
sdCursorAngle           = NaN(nBins,nSubj,nDays);
handAngle               = NaN(nBins,nSubj,nDays);
sdHandAngle             = NaN(nBins,nSubj,nDays);
reportAngle             = NaN(nBins,nSubj,nDays);
aimFixAngle             = NaN(nBins,nSubj,nDays);
implicitAngle           = NaN(nBins,nSubj,nDays);
implicitAngle_fromFix   = NaN(nBins,nSubj,nDays);
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
        targetZone = Results.fixAngles.day1.preview.targetZone;
    end
    
    %% Bin reaction time, cursor angle, hand angle, report angle, and fixation angle
    
    reportAngle_all = NaN(nBins,nDays);
    reportAngle_noOutliers = NaN(nBins,nDays);
    
    % loop over bins
    iBin = 1 : nTargets : nTrials+1;
    for b = 1 : length(iBin)-1
        % get all values in bin
        rt = Results.RT(iBin(b):iBin(b+1)-1,:);                 % reaction time
        ca = Results.hitAngle_cursor(iBin(b):iBin(b+1)-1,:);    % cursor hit angle
        ha = Results.hitAngle_hand(iBin(b):iBin(b+1)-1,:);      % hand hit angle
        ra_all = Results.explicitAngle(iBin(b):iBin(b+1)-1,:);  % reported aiming angle
        ra_no = Results.explicitAngle_outliersRemoved(iBin(b):iBin(b+1)-1,:);
        % get all aim fixation angles in bin
        fa = NaN;
        if rotBins(b)
            fa = Results.fixAngle_closestA(iBin(b):iBin(b+1)-1,:);
        elseif woBins(b)
            fa = Results.fixAngle_closestOppA(iBin(b):iBin(b+1)-1,:);
        end
        fa(:,sum(~isnan(fa))==1) = NaN;  % do not compute mean when there is only 1 value
        % calculate bin mean or median
        if meanOrMedian==1
            RT(b,s,:) = nanmean(rt);
            cursorAngle(b,s,:) = nanmean(ca);
            handAngle(b,s,:) = nanmean(ha);
            reportAngle_all(b,:) = nanmean(ra_all);
            reportAngle_noOutliers(b,:) = nanmean(ra_no);
            aimFixAngle(b,s,:) = nanmean(fa);
        elseif meanOrMedian==2
            RT(b,s,:) = nanmedian(rt);
            cursorAngle(b,s,:) = nanmedian(ca);
            handAngle(b,s,:) = nanmedian(ha);
            reportAngle_all(b,s,:) = nanmedian(ra_all);
            reportAngle_noOutliers(b,s,:) = nanmedian(ra_no);
            aimFixAngle(b,s,:) = nanmedian(fa);
        end
        % calculate bin standard deviation
        sdRT(b,s,:) = nanstd(rt);
        sdCursorAngle(b,s,:) = nanstd(ca);
        sdHandAngle(b,s,:) = nanstd(ha);
    end
    iBin = iBin(1:end-1)';
    
    % create separate aim fixation variables for rotation and washout?
    
    %% Check removed outliers in reported aiming angles
    
    keepOutliers = false;
    outlierDiff = ~isnan(reportAngle_noOutliers) & ...
        (reportAngle_all ~= reportAngle_noOutliers);
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
        reportAngle(:,s,:) = reportAngle_all;
        keyboard
    else
        reportAngle(:,s,:) = reportAngle_noOutliers;
    end
    
    %% Calculate implicit angle
    
    implicitAngle(:,s,:) = handAngle(:,s,:)-reportAngle(:,s,:);
    implicitAngle_fromFix(:,s,:) = handAngle(:,s,:)-aimFixAngle(:,s,:);
    
    %% Trials and bins in which an aiming strategy was present
    
    % first rotation bin in which an aiming strategy was present
    firstRotBin = find(rotBins,1);
    firstAimFixBin(s) = nanmax([find(aimFixAngle(firstRotBin:end,s,1)<targetZone(1),1) NaN]);
    firstAimReportBin(s) = nanmax([find(reportAngle(firstRotBin:end,s,1)<targetZone(1),1) NaN]);
    
    %% Plots
    
    if createPlots
        %
    end
    
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
