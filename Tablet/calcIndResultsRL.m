function calcIndResultsRL(projectPath,processData,createPlots,savePlots)
% calcIndResultsRL Calculate individual participant results in RL
% experiments.
%
% calcIndResultsRL(projectPath,processData,createPlots,savePlots) 
% loads the individual data in projectPath, and processes and saves the 
% individual data when processData=TRUE. Individual results are plotted 
% when createPlots=TRUE, and plots are saved when savePlots=TRUE. If the 
% data has been processed before, plots can be created using the saved data
% (processData=FALSE).

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

close all

% set defaults
if nargin==0
    projectPath = [];
    processData = true;
    createPlots = true;
    savePlots = false;
end

% select project data path
if isempty(projectPath)
    projectPath = '/Users/anouk/Documents/ExpsTablet/';
    expFolder = selectFiles([projectPath '*RL*'],'folders');
    projectPath = [projectPath expFolder.name '/'];
end
% check if we are at the right level
while ~exist([projectPath '/1_RawData/'],'dir');
    expFolder = selectFiles(projectPath,'folders');
    projectPath = [projectPath expFolder.name '/'];
end
expName = getFolderName(projectPath); % experiment name
% define path for input and output data
dataPath = [projectPath '/2_ProcessedData/'];
saveToPath = [projectPath '/3_Results/'];
if ~exist(saveToPath,'dir')
    mkdir(saveToPath)
end
saveFigsToPath = [saveToPath 'Figures/'];
if ~exist(saveFigsToPath,'dir')
    mkdir(saveFigsToPath)
end

% select subject data files
if processData
    cd(dataPath);
else
    cd(saveToPath);
end
subjFiles = selectFiles('S*.mat');
nSubj = length(subjFiles);

%% Loop over subjects

for s = 1 : nSubj
    fprintf('\n')
    
    if processData
        
        % get data file
        load([dataPath subjFiles(s).name]);
        disp(['Loaded ' subjFiles(s).name])
        
        % general info
        blockNo = Exp.blockNo;
        nBlocks = max(blockNo);
        nTrials_block = zeros(nBlocks,1);
        for b = 1 : nBlocks
            currBlock = blockNo==b';
            nTrials_block(b) = sum(currBlock);
        end
        
        % get relevant variables
        score = D.score;
        pathChange = D.pathChange;
        RT = D.iTargetGoLeaveCompleteEnd(:,3)-D.iTargetGoLeaveCompleteEnd(:,2);
        MT = D.iTargetGoLeaveCompleteEnd(:,4)-D.iTargetGoLeaveCompleteEnd(:,3);
        
        % get relation between score on trial t and path change on trial t+1
        score_ = score(blockNo>1); % remove practice block
        pathChange_ = pathChange(blockNo>1);
        score_pathChange = [score_(1:end-1) pathChange_(2:end)];
        score_pathChange = score_pathChange(~isnan(score_pathChange(:,2)),:); % remove last/first trials of block
        
        %% Loop over blocks to get scores, path changes, RT and MT
        
        % sort blocks based on number of curves and curvature
        pathDir = zeros(size(Exp.targetPathCurvature));
        pathDir(Exp.targetPathCurvature<0) = -1; % leftward curves
        pathDir(Exp.targetPathCurvature>0) = 1; % rightward curves
        expSpecs = [Exp.targetPathNumberOfCurves abs(Exp.targetPathCurvature) pathDir];
        [~,blockOrder] = sortrows(expSpecs(2:end,:));
        blockOrder = [1; blockOrder+1]; % start with baseline block
        % to do: save original block order to test for learning across blocks
        
        % preallocate
        allScores = NaN(nBlocks,max(nTrials_block));
        allPathChanges = NaN(nBlocks,max(nTrials_block));
        allRT = NaN(nBlocks,max(nTrials_block));
        allMT = NaN(nBlocks,max(nTrials_block));
        
        % create matrix with scores per block
        for b = 1 : nBlocks
            currBlock = blockOrder(b)==blockNo;
            allScores(b,1:sum(currBlock)) = score(currBlock);
            allPathChanges(b,1:sum(currBlock)) = pathChange(currBlock);
            allRT(b,1:sum(currBlock)) = RT(currBlock);
            allMT(b,1:sum(currBlock)) = MT(currBlock);
        end
        
        %% Save individual results
        
        ExpDetails = Exp;
        ExpDetails.nTrials_block = nTrials_block;
        
        Data = D;
        
        Results.blockNo = blockOrder;
        Results.targetPathCurvature = abs(Exp.targetPathCurvature(blockOrder));
        Results.targetPathDirection = pathDir(blockOrder);
        Results.targetPathNumberOfCurves = Exp.targetPathNumberOfCurves(blockOrder);
        Results.targetPathVisible = Exp.targetPathVisible(blockOrder);
        Results.score = allScores;
        Results.pathChange = allPathChanges;
        Results.RT = allRT;
        Results.MT = allMT;
        Results.score_pathChange = score_pathChange;
        
        % check if file exists
        fileName = [Exp.subjFolder '.mat'];
        if exist([saveToPath fileName],'file') == 2
            disp(['A file named ' fileName ' already exists.'])
            overwrite = input('Do you want to overwrite it? Yes(1) or no(0): ');
        else
            overwrite = 1;
        end
        
        % save file
        if overwrite == 1
            save([saveToPath fileName],'ExpDetails','Data','Results');
            disp(['Saved ' saveToPath fileName])
        else
            disp('Results not been saved')
        end
        
    elseif ~processData % load data
        load([saveToPath subjFiles(s).name])
        disp(['Loaded ' subjFiles(s).name])
    end
    
    %% Plots
    
    if createPlots
        
        %% Details for plots
        
        % blocks and trials
        nBlocks = max(Results.blockNo);
        nTrials_block = ExpDetails.nTrials_block;
        
        % specify subplots, block 1 = baseline
        nRow = 5;
        nCol = ceil(nBlocks/nRow);
        % create subblocks for trajectory plot if more than 50 trials per block
        if any(nTrials_block(2:end)>50)
            nPlots_sub = ceil(sum(nTrials_block(2:end))/50)+1;
            nCol_sub = ceil(nPlots_sub/nRow_sub);
        else
            nPlots_sub = nBlocks;
            nCol_sub = nCol;
        end
        
        %% Plot scores per block
        
        % create subplot with scores per block
        fig1 = scaledFigure(nCol*.6,nRow*.5);
        colors = get(gca,'colororder');
        for b = 1 : nBlocks
            currBlock = ExpDetails.blockNo==b;
            currScores = Data.score(currBlock,:);
            subplot(nRow,nCol_sub,b); hold on
            p = plot(currScores(:,1),'.-','color',colors(ExpDetails.targetPathNumberOfCurves(b),:,:));
            ylim([0 100]); xlabel('Trial'); ylabel('Score');
            str = sprintf('Block %d: path curvature = %1.2f',b,ExpDetails.targetPathCurvature(b));
            title(str)
            if b==1 && max(ExpDetails.targetPathNumberOfCurves)>1
                text(1,60,'Single curve','color',colors(1,:)) % legend
                text(1,40,'Double curve','color',colors(2,:))
            end
        end
        suplabel(['Reinforcement learning raw scores - ' ExpDetails.subjFolder],'t'); % grand title
        
        if savePlots
            saveFigAsPDF([saveFigsToPath 'RL_rawScores_' ExpDetails.subjFolder],10)
        end
        
        %% Plot cursor trajectories
        
        % plot data in xy coordinates
        fig2 = scaledFigure(nCol_sub*.6,nRow*.5);
        nTrials_baseline = nTrials_block(1);
        for i = 1 : nPlots_sub
            if i==1 % baseline
                iUseTrials = find(ExpDetails.blockNo==1);
                n = length(iUseTrials);
                colors = createColorGradient([0 0 1],[0 1 0],n); % blue to green
            elseif nPlots_sub==nBlocks % 1 subplot per block
                iUseTrials = find(ExpDetails.blockNo==i);
                n = length(iUseTrials);
                colorgradient = createColorGradient([0 0 1],[0 1 0],n);
                colors(iUseTrials,:) = colorgradient;
                colormap(colors)
            else % 50 trial intervals
                iUseTrials = ((i-2)*50+1 : (i-1)*50)' + nTrials_baseline;
                colorgradient = createColorGradient([0 0 1],[0 1 0],n);
                colors = [colors; colorgradient];
            end
            b = ExpDetails.blockNo(iUseTrials(1));
            
            % plot
            subplot(nRow,nCol_sub,i)
            hold on
            % plot start, target path, and target
            ps = plot(0,0,'ko');
            pv = plot(ExpDetails.xVisiblePath_plot(:,b),ExpDetails.yVisiblePath_plot(:,b),...
                'k--','linewidth',2);
            pp = plot(ExpDetails.xTargetPath_plot(:,b),ExpDetails.yTargetPath_plot(:,b),...
                'k','linewidth',2);
            pt = rectangle('Position',[ExpDetails.stim.targetXY-0.5*ExpDetails.stim.targetWidthHeight,...
                ExpDetails.stim.targetWidthHeight]);
            
            % plot cursor trajectories
            pc_d = [];
            for t = iUseTrials' % must be row vector to loop
                xyCursor = Data.xyCursor{t};
                int = Data.iTargetGoLeaveCompleteEnd(t,1):Data.iTargetGoLeaveCompleteEnd(t,4);
                if ~isnan(int)
                    if t==iUseTrials(1) || t==iUseTrials(end)
                        % use dashed line for first and last trial
                        pc_d(t) = plot(xyCursor(int,1),xyCursor(int,2),'--','color',colors(t,:));
                    else
                        pc = plot(xyCursor(int,1),xyCursor(int,2),'-','color',colors(t,:));
                    end
                end
            end
            
            % axes
            axis equal
            targetXY = ExpDetails.stim.targetXY;
            axis([-targetXY(2) targetXY(2) -10 targetXY(2)+10])
            if i==1
                legend([pp pc_d(1) pc_d(end)],...
                    {'target path','cursor - first trial','cursor - last trial'},'location','west');
                legend('boxoff')
            end
            xlabel('X (mm)'); ylabel('Y (mm)')
            str = sprintf('Block %d trials %d to %d',b,...
                ExpDetails.trialNo(iUseTrials(1)),ExpDetails.trialNo(iUseTrials(end)));
            title(str);
            hold off
        end
        sup = suplabel(['Cursor trajectories - ' ExpDetails.subjFolder],'t'); % grand title
        
        if savePlots
            saveFigAsPDF([saveFigsToPath 'RL_trajectories_' ExpDetails.subjFolder],10)
        end
        
        %% Plot relation between score and path change on next trial
        
        fig3 = figure; movegui(fig3,'northeast');
        scatterWithLinearFit(Results.score_pathChange(:,1),Results.score_pathChange(:,2),...
            [0 100 0 max(pathChange)+5]);
        set(gca,'xtick',0:20:100)
        set(gca,'ytick',0:20:max(pathChange)+5)
        xlabel('Score on trial t')
        ylabel('Path change on trial t+1')
        title(['Relation between score and path change - ' ExpDetails.subjFolder],...
            'interpreter','none')
        
        if savePlots
            saveFigAsPDF([saveFigsToPath 'RL_correlation_score_pathChange_' ExpDetails.subjFolder],12)
        end
        
        %% End
        
        disp('Press continue to go to next subject.')
        keyboard
        
        % close figures
        close([fig1 fig2 fig3])
        
    end % if createPlots
    
end % loop over subjects

%% Save copy of code used to process data

mFilePath = mfilename('fullpath');
saveCopyOfCode(mFilePath,saveToPath)
