function plotIndResultsVMRgaze(Results,ExpDetails,savePlots,saveToPath)
% plotIndResultsVMRgaze  Plot individual participant results calculated in calcIndResultsVMRgaze.m

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

%% Details for plots

% days
nDays = size(ExpDetails,2);
days = {'day1','day2'};
% subject
if nDays == 1
    subjName = ExpDetails(1).subjFolder;
else
    subjName_ = ExpDetails(1).subjName;
    subjName = ExpDetails(1).subjFolder(1:end-5);
    if ~strcmp(subjName_,subjName) % temp
        keyboard
    end
end
% blocks, trials and bins
nBlocks = max(Results.blockNo);
nTrials = sum(~isnan(Results.trialNo(:,1:nDays)));
maxnTrials = max(nTrials);
if isfield(ExpDetails(1),'breakTrials')
    breakTrials = [ExpDetails.breakTrials];
else
    breakTrials = [];
end
iBin = Results.trialNo_bins(:,end)+3;
rotBins = false(length(iBin),nDays);
% sampling
fs = ExpDetails(1).setup.fs;
tsample = 1000/fs;
analyzedGazeData = any(Results.analyzedGazeData);
plotGaze = true;
% specs per day
cursorRotation = Results.cursorRotation;
rotationOn = cursorRotation~=0 & ~isnan(cursorRotation);
iNewBlock = NaN(max(nBlocks),nDays);
nTrials_blocks = NaN(max(nBlocks),nDays);
previewDur = NaN(maxnTrials,nDays);
for d = 1 : nDays
    % number of trials
    n = length(ExpDetails(d).trialNo);
    % rotation
    rotBins(:,d) = cursorRotation(iBin,d)~=0;
    dRot = find(diff([0;rotationOn(:,d);0]));
    iRotationOnOff(1:length(dRot),d) = dRot-0.5;
    % blocks
    for b = 1 : nBlocks(d);
        iNewBlock(b,d) = nanmax([find(Results.blockNo(:,d)==b,1) NaN]);
        nTrials_blocks(b,d) = sum(Results.blockNo(:,d)==b);
    end
    % target preview
    previewDur(1:n,d) = ExpDetails(d).timing.previewDur;
end
nRot = length(iRotationOnOff)/2;
vmr = cursorRotation(find(rotationOn,1));

% reporting
trialType = Results.trialType;
reportTrials = strcmp(trialType,'report');
plotReport = false;
if any(reportTrials(:))
    plotReport = true;
end
% gaze
if analyzedGazeData
    pearlAngles = mean(Results.percTimeFix.day1.pearlBinAngles);
else
    plotGaze = false;
end
% RT
plotRT = false;
if any(previewDur==0)
    plotRT = true;
end

%% Open figures

%close all
fig1 = scaledFigure(0.5+1.0*plotGaze+0.5*plotReport+0.5*plotRT,0.8*nDays);  % raw angles [hand explicit fix-preview fix-reach]
fig2 = scaledFigure(0.5+0.5*plotGaze+1*plotReport+0.5*plotRT,0.8*nDays);    % binned angles - subplots
fig2b = scaledFigure(1,0.8*nDays);                                          % binned angles - overlayed
if plotGaze
    fig3 = scaledFigure(1,0.8*nDays);                                       % distribution of fixation angles per block
    if ~isempty(Results.probFix_blocks)
        nSubblocks = length(fieldnames(Results.probFix_subblocks.day1));
    else
        nSubblocks = 1;
    end
    fig4 = scaledFigure(2,0.8*nDays);               % probability of fixation per block
    fig4b = scaledFigure(2/3*nSubblocks,0.8*nDays); % probability of fixation per rotation subblock
    fig5 = scaledFigure(1,0.8*nDays);               % histogram of fixation angles for sequential fixations
end
colors = get(gca,'colororder');
fadedColors = colors+(1-colors)*0.6;

%% Plot raw hand angles, reported aim angles and fixation angles

dataToPlot = {Results.hitAngle_hand};
figTitles = {'Learning'};
yLabels = {'Endpoint hand angle (deg)'};
if plotReport
    if ~isfield(Results,'explicitAngle_outliersRemoved')
        Results.explicitAngle_outliersRemoved = Results.explicitAngle;
    end
    dataToPlot = [dataToPlot {Results.explicitAngle_outliersRemoved}];
    figTitles = [figTitles {'Explicit learning'}];
    yLabels = [yLabels {'Reported aim angle (deg)'}];
end
if plotGaze
    dataToPlot = [dataToPlot {cat(3,Results.fixAngle_closestA,Results.fixAngle_closestOppA),NaN(maxnTrials,2)}];
    figTitles = [figTitles {'Fixation during target preview','Fixation during reach'}];
    yLabels = [yLabels {'Fixation angle (deg)','Fixation angle (deg)'}];
end
if plotRT
    dataToPlot = [dataToPlot {Results.RT*1000}]; % in ms
    figTitles = [figTitles {'Reaction time'}];
    yLabels = [yLabels {'Reaction time (ms)'}];
end

% plot raw hand and explicit angles, plot all fixations,
% with fixations closest to aimpoint in a different color
figure(fig1); clf
nCol = length(dataToPlot);
for d = 1 : nDays
    for c = 1 : nCol
        % create subplot with correct size
        h(c) = subplot(nDays,nCol,d*nCol-nCol+c); pb = pbaspect; hold on
        pb = pbaspect; pbaspect([pb(1) 0.95*pb(2) pb(3)]); h(c).Position(2) = h(c).Position(2)-0.025*pb(2);
        % set axes
        xlim([0 maxnTrials+1]); set(gca,'XTick',0:80:maxnTrials)
        if ~isempty(strfind(yLabels{c},'deg')); ylim([-95 50]); set(gca,'YTick',-90:45:45)
        elseif ~isempty(strfind(yLabels{c},'ms')); ylim([0 max(ExpDetails(1).timing.maxRT)*1000]); end
        yl = ylim;
        % color background when rotation is on, works in R2015b but not in 2017a
        for r = 1 : nRot
            x = iRotationOnOff(r*2-1:r*2,d);
            a = area([x; flipud(x)],[yl(1) yl(1) yl(2) yl(2)],'LineStyle','none');
            a.FaceColor = [0.95 0.95 0.95];
            plot(x,-[vmr vmr],'k') % draw line at hand target angle
        end
        plot([0 maxnTrials],[0 0],'k'); % draw line at 0
        % plot all fixation angles
        if strcmp(figTitles(c),'Fixation during target preview') && plotGaze
            plot(Results.fixAngles.(days{d}).preview.all,'.','color',fadedColors(d,:));
            text(70,-72,'all','color',fadedColors(d,:))
            text(70,-80,['closest to +/-' num2str(vmr)],'color',colors(d,:))
        elseif strcmp(figTitles(c),'Fixation during reach') && plotGaze
            plot(Results.fixAngles.(days{d}).reach.all,'.','color',fadedColors(d,:));
            text(70,-72,'all','color',fadedColors(d,:))
        end
        % plot data
        plot(squeeze(dataToPlot{c}(:,d,:)),'.','color',colors(d,:));
        vertline(breakTrials(:,d),'g:')
        vertline([iNewBlock(:,d)-0.5; nTrials(d)+0.5],'k:')
        hold off
        % label axes and add title
        xlabel('Trial')
        ylabel(yLabels{c})
        if d==1
            title(figTitles{c})
        end
    end
end
suplabel(['Raw angles - ' subjName],'t'); % grand title

% save
if savePlots
    saveFigAsPDF([saveToPath 'rawAngles_' subjName],12)
end

%% Plot binned angles - subplots

dataToPlot = {Results.mnBinnedHitAngle_hand};
sdDataToPlot = {Results.sdBinnedHitAngle_hand};
figTitles = {'Learning'};
yLabels = {'Endpoint hand angle (deg)'};
legends = {[]};
if plotReport
    if ~isfield(Results,'mnBinnedExplicitAngle_outliersRemoved')
        Results.mnBinnedExplicitAngle_outliersRemoved = Results.binnedExplicitAngle;
    end
    Results.mnBinnedImplicitAngle_fix(~rotBins) = NaN;
    dataToPlot = [dataToPlot {Results.mnBinnedExplicitAngle_outliersRemoved,...
        cat(3,Results.mnBinnedImplicitAngle,Results.mnBinnedImplicitAngle_fix)}];
    sdDataToPlot = [sdDataToPlot {[],cat(3,[],[])}];
    figTitles = [figTitles {'Explicit learning','Implicit learning'}];
    yLabels = [yLabels {'Reported aim angle (deg)','Hand minus reported aim angle (deg)'}];
    legends = [legends {[],{'From report'}}];
    %legends = [legends {[],{'From report','From fixation'}}];
end
if plotGaze
    Results.mnBinnedFixAngle_closestA(~rotBins) = NaN;
    dataToPlot = [dataToPlot {cat(3,Results.mnBinnedFixAngle_closestA,...
        Results.mnBinnedFixAngle_closestOppA)}];
    sdDataToPlot = [sdDataToPlot {cat(3,Results.sdBinnedFixAngle_closestA,...
        Results.sdBinnedFixAngle_closestOppA)}];
    figTitles = [figTitles {'Preview fixation closest to hand target'}];
    yLabels = [yLabels {'Fixation angle (deg)'}];
    legends{end} = {'From report','From fixation'};
    legends = [legends {{['Closest to -' num2str(vmr)],['Closest to +' num2str(vmr)]}}];
end
if plotRT
    dataToPlot = [dataToPlot {Results.mnBinnedRT*1000}]; % in ms
    sdDataToPlot = [sdDataToPlot {Results.sdBinnedRT*1000}];
    figTitles = [figTitles {'Reaction time'}];
    yLabels = [yLabels {'Reaction time (ms)'}];
    legends = [legends {[]}];
end
colors = cat(3,colors,fadedColors); % concatenate for looping

figure(fig2); clf;
nCol = length(dataToPlot);
for d = 1 : nDays
    for c = 1 : nCol
        % create subplot with correct size
        h(c) = subplot(nDays,nCol,d*nCol-nCol+c); hold on
        pb = pbaspect; pbaspect([pb(1) 0.95*pb(2) pb(3)]);
        h(c).Position(2) = h(c).Position(2)-0.025*pb(2);
        % set axes
        xlim([0 maxnTrials+1]); set(gca,'XTick',0:80:maxnTrials)
        if ~isempty(strfind(yLabels{c},'deg'));
            ylim([-95 50]); set(gca,'YTick',-90:45:45)
        elseif ~isempty(strfind(yLabels{c},'ms'));
            ylim([0 max(ExpDetails(1).timing.maxRT)*1000]);
        end
        yl = ylim;
        % color background when rotation is on
        for r = 1 : nRot
            x = iRotationOnOff(r*2-1:r*2,d);
            a = area([x; flipud(x)],[yl(1) yl(1) yl(2) yl(2)],'LineStyle','none');
            a.FaceColor = [0.95 0.95 0.95];
            plot(x,-[vmr vmr],'k') % draw line at hand target angle
        end
        plot([0 maxnTrials],[0 0],'k'); % draw line at 0
        % plot data
        for i = size(dataToPlot{c},3):-1:1
            p(i) = plot(iBin,dataToPlot{c}(:,d,i),'o','color',colors(d,:,i));
            if ~isempty(sdDataToPlot{c}(:,:,i))
                errorb(iBin,dataToPlot{c}(:,d,i),sdDataToPlot{c}(:,d,i),'barwidth',0,'color',colors(d,:,i))
            end
        end
        hold off
        % add labels, title, and legend
        vertline(breakTrials(:,d),'g:')
        vertline([iNewBlock(:,d)-0.5; nTrials(d)+0.5],'k:')
        xlabel('Trial')
        ylabel(yLabels{c})
        if d==1
            title(figTitles{c})
            if ~isempty(legends{c})
                legend([p(1) p(2)],legends{c},'location','southwest'); legend('boxoff')
            end
        end
    end
end
suplabel(['Binned angles - ' subjName],'t');

% save
if savePlots
    saveFigAsPDF([saveToPath 'binnedAngles_' subjName],12)
end

colors = colors(:,:,1); % reset

%% Plot binned angles - single plot overlayed

figure(fig2b); clf
for d = 1 : nDays
    % create subplot with correct size
    h(c) = subplot(nDays,1,d); hold on
    pb = pbaspect; pbaspect([pb(1) 0.95*pb(2) pb(3)]);
    h(c).Position(2) = h(c).Position(2)-0.025*pb(2);
    % set axes
    xlim([0 maxnTrials]); set(gca,'XTick',0:80:maxnTrials)
    ylim([-60 20]); yl = ylim; set(gca,'Ytick',-60:15:45)
    % color background when rotation is on
    for r = 1 : nRot
        x = iRotationOnOff(r*2-1:r*2,d);
        a = area([x; flipud(x)],[yl(1) yl(1) yl(2) yl(2)],'LineStyle','none');
        a.FaceColor = [0.95 0.95 0.95];
        plot(x,-[vmr vmr],'k') % draw line at hand target angle
    end
    plot([0 maxnTrials],[0 0],'k'); % draw line at 0
    % plot data
    pi = plot(iBin,Results.mnBinnedImplicitAngle(:,d),'o-','color',colors(5,:));
    ph = plot(iBin,Results.mnBinnedHitAngle_hand(:,d),'o-','color',colors(1,:));
    pe = plot(iBin,Results.mnBinnedExplicitAngle_outliersRemoved(:,d),'o-','color',colors(2,:));
    %pif = plot(iBin,Results.mnBinnedImplicitAngle_fix(:,d),'o-','color',fadedColors(5,:));
    pf = plot(iBin,Results.mnBinnedFixAngle_closestA(:,d),'o-','color',colors(4,:));
    hold off
    % add lables and title
    vertline(breakTrials(:,d),'g:')
    vertline([iNewBlock(:,d)-0.5; nTrials(d)+0.5],'k:')
    xlabel('Trial')
    ylabel('Angle (deg)')
    if d==1
        title(['Binned angles - ' subjName],'interpreter','none')
        legend([ph pe pi pf],{'Hand','Explicit','Implicit','Fixation'},'location','east');
        legend('boxoff')
    end
end

if savePlots
    saveFigAsPDF([saveToPath 'binnedOverlayedAngles_' subjName],12)
end

%% Plot distributions of fixation angles per block

if any(Results.analyzedGazeData)
    
    figure(fig3); clf
    % plot gaze angle relative to target - target preview and reach
    for d = 1 : nDays
        b = 1 : size(Results.percTimeFix_blocks.(days{d}).preview.atPearl,1);
        % skip baseline-report block and short blocks
        iBRblock = find(strncmpi(Results.percTimeFix_blocks.(days{d}).blockNames,'Baseline+report',15));
        iShortBlock = find(nTrials_blocks(:,d)<=16);
        if ~isempty([iBRblock iShortBlock]); b([iBRblock iShortBlock])=[]; end ;
        nCol = length(b);
        for c = 1 : nCol
            h(c) = subplot(nDays,nCol,d*nCol-nCol+c); hold on
            pb = pbaspect; pbaspect([pb(1) 0.95*pb(2) pb(3)]); h(c).Position(2) = h(c).Position(2)-0.025*pb(2);
            plot(pearlAngles,Results.percTimeFix_blocks.(days{d}).preview.atPearl(b(c),:),...
                '.-','color',colors(d,:)); % preview
            plot(pearlAngles,Results.percTimeFix_blocks.(days{d}).reach.atPearl(b(c),:),...
                '.-','color',fadedColors(d,:)); % reach
            axis([-90 90 0 80]); set(gca,'Xtick',-90:45:90)
            vertline(-90:22.5:90)
            title(Results.percTimeFix_blocks.(days{d}).blockNames(b(c)))
            if c == 1
                text(-80,75,'Preview','Color',colors(d,:))
                text(-80,70,'Reach','Color',fadedColors(d,:))
            end
        end
    end
    supx = suplabel('Angle relative to target (deg)','x');
    supx.Position = [supx.Position(1) supx.Position(2)+0.04 supx.Position(3:4)];
    supy = suplabel('% Fixation time','y');
    supy.Position = [supy.Position(1)+0.03 supy.Position(2:4)];
    suplabel(['Fixation angles during preview and reach - ' subjName],'t');
    clear supx supy
    
    if savePlots
        saveFigAsPDF([saveToPath 'fixAngles_' subjName],12)
    end
    
    %% Plot probability of fixation in start, target and aimpoint area - blocks
    
    if ~isempty(Results.probFix_blocks)
        
        % define zones
        targetBin = find(pearlAngles==0);
        targetBins = targetBin-1:targetBin+1;
        aimZoneBins = find(pearlAngles>-vmr,1):(targetBins(1)-1);
        labels = {'Start area','Target area',...
            ['Area between -' num2str(vmr) ' and target']};
        
        %
        probFix = Results.probFix_blocks;
        tPhase = probFix.intervals*tsample;
        blocks = fieldnames(probFix.day1);
        b = 1 : length(blocks);
        b(iShortBlock) = []; % skip baseline-report block
        nCol = length(b);
        
        % plot
        figure(fig4); clf
        for d = 1 : nDays
            for c = 1 : nCol
                % compute
                probFixStart = probFix.(days{d}).(blocks{b(c)}).start;
                probFixPearls_target = sum(probFix.(days{d}).(blocks{b(c)}).pearls(:,targetBins),2);
                probFixPearls_aim = sum(probFix.(days{d}).(blocks{b(c)}).pearls(:,aimZoneBins),2);
                t = 0:tsample:(length(probFixStart)-1)*tsample;
                % plot
                h(c) = subplot(nDays,nCol,d*nCol-nCol+c); hold on
                pb = pbaspect; pbaspect([pb(1) 0.95*pb(2) pb(3)]);
                h(c).Position(2) = h(c).Position(2)-0.025*pb(2);
                ps=plot(t,probFixStart);
                pt=plot(t,probFixPearls_target);
                pa=plot(t,probFixPearls_aim);
                axis([0 tPhase(end)-20 0 1])
                vertline(tPhase,'k:')
                %xlabel('Normalized time (ms)')
                %ylabel('Probability of fixation');
                if d==1 && c==2
                    legend([ps pt pa],'Start area','Target area',...
                        'Area between -45 and target','location','northwest');
                    legend('boxoff')
                end
                if round(b(c))==b(c)
                    keyboard
                    %title(blockNames{b(c)})
                end
            end
        end
        supx = suplabel('Normalized time (ms)','x');
        supx.Position = [supx.Position(1) supx.Position(2)+0.04 supx.Position(3:4)];
        supy = suplabel('Probability of fixation','y');
        supy.Position = [supy.Position(1)+0.03 supy.Position(2:4)];
        suplabel(['Probability of fixating in start, target, and aimpoint area - ' subjName],'t');
        clear supx supy
        
        if savePlots
            saveFigAsPDF([saveToPath 'fixProbabilities_' subjName],12)
        end
        
    end
    
    %% Plot probability of fixation in start, target and aimpoint area -  rotation subblocks
    
    if isfield(Results,'probFix_subblocks') && ~isempty(Results.probFix_subblocks)
        
        % subblocks in rotation block
        probFix = Results.probFix_subblocks;
        tPhase = probFix.intervals*tsample;
        blocks = fieldnames(probFix.day1);
        b = 1 : length(blocks);
        nCol = length(b);
        
        % plot
        figure(fig4b); clf
        for d = 1 : nDays
            for c = 1 : nCol
                % compute
                probFixStart = probFix.(days{d}).(blocks{b(c)}).start;
                probFixPearls_target = sum(probFix.(days{d}).(blocks{b(c)}).pearls(:,targetBins),2);
                probFixPearls_aim = sum(probFix.(days{d}).(blocks{b(c)}).pearls(:,aimZoneBins),2);
                t = 0:tsample:(length(probFixStart)-1)*tsample;
                % plot
                h(c) = subplot(nDays,nCol,d*nCol-nCol+c); hold on
                pb = pbaspect; pbaspect([pb(1) 0.95*pb(2) pb(3)]);
                h(c).Position(2) = h(c).Position(2)-0.025*pb(2);
                ps=plot(t,probFixStart);
                pt=plot(t,probFixPearls_target);
                pa=plot(t,probFixPearls_aim);
                axis([0 tPhase(end)-20 0 1])
                vertline(tPhase,'k:')
                %xlabel('Normalized time (ms)')
                %ylabel('Probability of fixation');
                if d==1 && c==2
                    legend([ps pt pa],'Start area','Target area',...
                        'Area between -45 and target','location','northwest');
                    legend('boxoff')
                end
                title(blocks{c},'Interpreter','none')
            end
        end
        supx = suplabel('Normalized time (ms)','x');
        supx.Position = [supx.Position(1) supx.Position(2)+0.04 supx.Position(3:4)];
        supy = suplabel('Probability of fixation','y');
        supy.Position = [supy.Position(1)+0.03 supy.Position(2:4)];
        suplabel(['Probability of fixating in start, target, and aimpoint area - ' subjName],'t');
        clear supx supy
        
        if savePlots
            saveFigAsPDF([saveToPath 'fixProbabilitiesRot_' subjName],12)
        end
        
    end
    
    %% Histogram of fixation angles, for each sequential fixation
    
    if length(pearlAngles)>1
        pa = pearlAngles(25:41);
        
        figure(fig5); clf
        for d = 1 : nDays
            fixAngles_rot = Results.fixAngles.(days{d}).preview.all(rotationOn(:,d),:);
            fixAngles_rot = fixAngles_rot(~all(isnan(fixAngles_rot),2),:);
            % plot
            h = subplot(nDays,1,d);
            pb = pbaspect; pbaspect([pb(1) 0.95*pb(2) pb(3)]);
            h.Position(2) = h.Position(2)-0.025*pb(2);
            hist(fixAngles_rot(:,1:end-1),pa);
            % axes and legend
            set(gca,'Xtick',pa); set(gca,'XTicklabelrotation',45);
            pos = get(gca,'Position');
            pos = [pos(1) pos(2)+0.05*pos(4) pos(3) 0.95*pos(4)];
            set(gca, 'Position', pos)
            %xlim([min(pa) max(pa)])
            legend('fix1','fix2','fix3','fix4','fix5','fix6','fix7','fix8')
            legend('boxoff')
            if nDays>1
                title(days{d})
            end
        end
        supx = suplabel('Fixation angle','x');
        %supx.Position = [supx.Position(1) supx.Position(2)+0.04 supx.Position(3:4)];
        supy = suplabel('Frequency','y');
        supy.Position = [supy.Position(1)+0.03 supy.Position(2:4)];
        suplabel(['Histogram of first to last fixation angle during preview - ' subjName],'t');
        clear supx supy
        
        if savePlots
            saveFigAsPDF([saveToPath 'fixHistogram_' subjName],12)
        end
        
    end % if length(pearlAngles)>1
    
end % if any(Results.analyzedGazeData)

%% End

disp('Press continue to go to next subject.')
keyboard

% close figures
close([fig1 fig2 fig2b])
if plotGaze
    close([fig3 fig4 fig4b fig5])
end
