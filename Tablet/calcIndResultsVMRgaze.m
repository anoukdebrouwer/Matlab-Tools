function calcIndResultsVMRgaze(processData,createPlots,savePlots)
% calcIndResultsVMRgaze  Calculate individual participant results in VMR
% experiments with gaze tracking.

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

close all

% set defaults
if nargin==0
    processData = true;
    createPlots = true;
    savePlots = false;
end

% select experiment data folder
projectPath = '/Users/anouk/Documents/ExpsTablet/';
expFolder = selectFiles(projectPath,'folders');
while ~exist([projectPath expFolder.name '/1_RawData/'],'dir');
    projectPath = [projectPath expFolder.name '/'];
    expFolder = selectFiles(projectPath,'folders'); % look in subfolders
end
dataPath = [projectPath expFolder.name '/2_ProcessedData/'];
saveToPath = [projectPath expFolder.name '/3_Results/'];
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
if ~isempty(dir('*day*mat')) % check for 2-day experiments
    subjFiles = selectFiles('*day1*.mat');
else
    subjFiles = selectFiles('*.mat');
end
nSubj = length(subjFiles);

% load experiment details
detailsFile = selectFiles([projectPath expFolder.name '/ExpDetails*.mat'],'files');
load([projectPath expFolder.name '/' detailsFile.name])
nTrialBins = nTrials/8;
blocks = cellstr([repmat('block',nBlocks,1) num2str((1:nBlocks)')]);
days = {'day1','day2'};
if exist('breakTrials','var')
    if iscell(breakTrials)
        breakTrials = breakTrials{expNumber};
        if size(breakTrials,2)<nDays
            breakTrials = [breakTrials(:) breakTrials(:)];
        end
    end
else
    breakTrials = [];
end

if processData
    % open figure
    %fig0 = scaledFigure(1.5,1.5);
    
    % criteria for analysis
    if ~exist('maxMT','var')
        maxMT = 0.4;
    end
    minMaxPercDistance = [75 125]; % gaze distance range for valid fixations
    maxPercNoGaze = 50; % maximum % time with missing data per trial
    
    % length of intervals that data is resampled to in order to average
    % trial time courses across subjects
    if exist('resamplingIntervals','var')
        nsNew = resamplingIntervals; % [preview RT MT post] n samples
    else
        nsNew = [];
    end
end

%% Loop over subjects

for s = 1 : nSubj
    fprintf('\n')
    
    if processData
        
        % preallocate
        clear ExpDetails Data Results
        feedbackAndMTGood   = false(nTrials,nDays);
        noReport            = false(nTrials,nDays);
        analyzedGazeData    = false(nTrials,nDays);
        fixAngles           = [];
        percTimeFix         = [];
        percTimeFix_blocks  = [];
        probFix_blocks      = [];
        probFix_subblocks   = [];
        blockNumber         = NaN(nTrials,nDays);
        trialNumber         = NaN(nTrials,nDays);
        targetPreviewDur    = NaN(nTrials,nDays);
        RT                  = NaN(nTrials,nDays);
        MT                  = NaN(nTrials,nDays);
        cursorRotation      = NaN(nTrials,nDays);
        hitAngle_cursor     = NaN(nTrials,nDays);
        hitAngle_hand       = NaN(nTrials,nDays);
        explicitAngle               = NaN(nTrials,nDays);
        explicitAngle_outliersRemoved = NaN(nTrials,nDays);
        fixAngle_preview_closestA   = NaN(nTrials,nDays);
        fixAngle_preview_closestOppA = NaN(nTrials,nDays);
        iFixTarget_aimFix_rot       = NaN(nTrials,nDays);
        mnBinnedHitAngle_cursor     = NaN(nTrialBins,nDays);
        medBinnedHitAngle_cursor    = NaN(nTrialBins,nDays);
        sdBinnedHitAngle_cursor     = NaN(nTrialBins,nDays);
        mnBinnedHitAngle_hand       = NaN(nTrialBins,nDays);
        medBinnedHitAngle_hand      = NaN(nTrialBins,nDays);
        sdBinnedHitAngle_hand       = NaN(nTrialBins,nDays);
        mnBinnedRT                  = NaN(nTrialBins,nDays);
        medBinnedRT                 = NaN(nTrialBins,nDays);
        sdBinnedRT                  = NaN(nTrialBins,nDays);
        mnBinnedExplicitAngle       = NaN(nTrialBins,nDays);
        medBinnedExplicitAngle      = NaN(nTrialBins,nDays);
        sdBinnedExplicitAngle       = NaN(nTrialBins,nDays);
        mnBinnedExplicitAngle_outliersRemoved = NaN(nTrialBins,nDays);
        sdBinnedExplicitAngle_outliersRemoved = NaN(nTrialBins,nDays);
        mnBinnedFixAngle_preview_closestA = NaN(nTrialBins,nDays);
        medBinnedFixAngle_preview_closestA = NaN(nTrialBins,nDays);
        sdBinnedFixAngle_preview_closestA = NaN(nTrialBins,nDays);
        mnBinnedFixAngle_preview_closestOppA = NaN(nTrialBins,nDays);
        medBinnedFixAngle_preview_closestOppA = NaN(nTrialBins,nDays);
        sdBinnedFixAngle_preview_closestOppA = NaN(nTrialBins,nDays);
        %figure(fig0);clf
        
        %% Loop over days
        
        for d = 1 : nDays
            
            % get data file
            if nDays == 2
                fileName = [subjFiles(s).name(1:end-5) num2str(d) '.mat'];
            else
                fileName = subjFiles(s).name;
            end
            load([dataPath fileName]);
            fprintf('Loaded %s',fileName)
            Exp.breakTrials = breakTrials(:,d);
            ExpDetails(d) = Exp;
            Data(d) = D;
            nTrials = length(Exp.trialNo); n = nTrials;
            nBlocks = max(Exp.blockNo);
            
            % general info
            blockNumber(1:n,d) = Exp.blockNo;
            trialNumber(1:n,d) = Exp.trialNo;
            rotBlockNumber = blockNumber(find(Exp.cursorRotation~=0,1),d);
            washoutBlockNumber = find(strcmpi(blockNames,'washout'));
            subBlockNumber = blockNumber(:,d);
            nTrials_rot = sum(blockNumber(:,d)==rotBlockNumber);
            nPearlBins = length(pearlBinAngles);
            
            %% Remove error trials
            
            % trial type
            noReport(1:n,d) = ~ismember(Exp.trialType,'ReachDM_G') & ...
                ~ismember(Exp.trialType,'ReachDM_K') & ~ismember(Exp.trialType,'ReachDM_L');
            
            % find error trials
            MT_temp = D.tTargetGoLeaveRingEnd(:,4)-D.tTargetGoLeaveRingEnd(:,3);
            feedbackAndMTGood(1:n,d) = ismember(D.feedbackMessage,'Good') & MT_temp<0.4;
            
            % remove error trials from struct
            D_new = replaceErrorTrialsInStruct(D,~feedbackAndMTGood(1:n,d),false);
            D_new.explicit1Angle = D.explicit1Angle; % keep report values
            D_new.feedbackAndMTGood = feedbackAndMTGood(1:n,d);
            
            % display number of correct trials
            str = sprintf('%d early trials (online)',sum(strcmp(D.feedbackMessage,'Early')));
            disp(str)
            str = sprintf('%d late trials (online)',sum(strcmp(D.feedbackMessage,'Late')));
            disp(str)
            str = sprintf('%d trials with MT > %1.1f s ',sum(MT_temp>maxMT | strcmp(D.feedbackMessage,'Short')),maxMT);
            disp(str)
            str = sprintf('%d non-reporting trials with good timing',sum(noReport(1:n,d)&feedbackAndMTGood(1:n,d)));
            disp(str)
            str = sprintf('%d reporting trials with good timing',sum(~noReport(1:n,d)&feedbackAndMTGood(1:n,d)));
            disp(str)
            
            %% Analyze gaze data
            
            noGazeData = cell2mat(cellfun(@(x) all(isnan(x(:))),D_new.xyGaze,'UniformOutput',false));
            %noGazeData = true;
            if ~all(noGazeData)
                [~,ind] = max(abs(ExpDetails(1).cursorRotation));
                vmr = ExpDetails(1).cursorRotation(ind);
                [percTimeFix_temp,fixAngles_temp,D_new] = analyzeVMRgazeData(Exp,D_new,minMaxPercDistance,nsNew,true);
                analyzedGazeData(1:n,d) = D_new.analyzedGazeData;
                
                % fixation angle closest to target and aimpoint
                iRows = (1:nTrials)';
                for st = 1 : length(states)
                    targetZone = fixAngles_temp.targetZone;
                    fixAngles_temp.(states{st}).closestT = NaN(nTrials,1);
                    fixAngles_temp.(states{st}).closestA = NaN(nTrials,1);
                    fixAngles_temp.(states{st}).closestOppA = NaN(nTrials,1);
                    fixAngles_temp.(states{st}).nCloseA_block = NaN(max(Exp.blockNo),1);
                    % remove columns without fixations
                    anyFix = any(~isnan(fixAngles_temp.(states{st}).all));
                    fixAngles_temp.(states{st}).all = fixAngles_temp.(states{st}).all(:,anyFix);
                    fixAngles_temp.(states{st}).iOnset_all = fixAngles_temp.(states{st}).iOnset_all(:,anyFix);
                    a = fixAngles_temp.(states{st}).all;
                    iOn = fixAngles_temp.(states{st}).iOnset_all;
                    % closest to target
                    [~,iMin] = nanmin(abs(a),[],2);
                    iMin = sub2ind(size(a),iRows,iMin); % get linear index
                    aClosest = a(iMin);
                    inTargetZone = aClosest>targetZone(1) & aClosest<targetZone(2);
                    fixAngles_temp.(states{st}).closestT(inTargetZone) = aClosest(inTargetZone);
                    % closest to aimpoint
                    [~,iMin] = nanmin(abs(a - (-vmr)),[],2);
                    iMin = sub2ind(size(a),iRows,iMin); % get linear index
                    aClosest = a(iMin);
                    iOnClosest = iOn(iMin);
                    inAimZone = aClosest>-90 & aClosest<targetZone(1);
                    fixAngles_temp.(states{st}).closestA(inAimZone) = aClosest(inAimZone);
                    fixAngles_temp.(states{st}).iOnset_closestA(inAimZone) = iOnClosest(inAimZone);
                    for b = 1 : max(Exp.blockNo)
                        currBlock = b==Exp.blockNo;
                        closestA_currBlock = fixAngles_temp.(states{st}).closestA(currBlock);
                        fixAngles_temp.(states{st}).nCloseA_block(b,1) = sum(~isnan(closestA_currBlock));
                        fixAngles_temp.nAnalyzed_block(b,1) = sum(D_new.analyzedGazeData(currBlock));
                    end
                    % closest to opposite aimpoint
                    [~,iMin] = nanmin(abs(a - vmr),[],2);
                    iMin = sub2ind(size(a),iRows,iMin); % get linear index
                    aClosest = a(iMin);
                    inAimZone = aClosest>targetZone(2) & aClosest<90;
                    fixAngles_temp.(states{st}).closestOppA(inAimZone) = aClosest(inAimZone);
                end
                
                % create struct with separate fields for each day
                fixAngles.(days{d}) = fixAngles_temp;
                percTimeFix.(days{d}) = percTimeFix_temp;
                
            end % if ~all(noGazeData)
            
            %% Average percentage of time that fixation is in pearl bin within blocks
            
            % create subblocks of 160 trials within rotation block
            if nTrials_rot>=160
                nSubBlocks_rot = nTrials_rot/160;
                for bsub = 1 : nSubBlocks_rot
                    subBlock = blockNumber(:,d)==rotBlockNumber & trialNumber(:,d)>nTrials_rot/nSubBlocks_rot*(bsub-1) &...
                        trialNumber(:,d)<=nTrials_rot/nSubBlocks_rot*bsub;
                    uniSubBlock_rot(bsub) = subBlockNumber(find(subBlock,1))+0.1*bsub;
                    subBlockNumber(subBlock) = uniSubBlock_rot(bsub);
                end
                subBlocks_rot = cellstr([repmat('subblock',length(uniSubBlock_rot),1) num2str(uniSubBlock_rot(:))]);
                subBlocks_rot = strrep(subBlocks_rot,'.','_');
            end
            
            if ~all(noGazeData)
                
                % for each block: average percentage of time that gaze direction is at pearl bin
                percTimeFix_blocks.(days{d}).blockNames = ExpDetails(d).blockNames;
                for st = 1 : length(states)
                    for b = 1 : nBlocks
                        currBlock = b == blockNumber(:,d);
                        percTimeFix_blocks.(days{d}).(states{st}).atPearl(b,:) = ...
                            nanmean(percTimeFix.(days{d}).(states{st}).atPearl(currBlock,:,1));
                        percTimeFix_blocks.(days{d}).(states{st}).atStart(b,:) = ...
                            nanmean(percTimeFix.(days{d}).(states{st}).atStart(currBlock));
                    end
                end
                
                %% Compute probability of fixating the target or aimpoint zone over time
                % for experiments with target preview only
                
                if ~isempty(nsNew) && ~isfield(D_new,'dirGaze_meanFix_resampled')
                    
                    % concatenate resampled intervals
                    dirGaze_resampled = cat(1,D_new.dirGaze_meanFix_resampled{:});
                    distGaze_resampled = cat(1,D_new.distGaze_meanFix_resampled{:});
                    nsTotal = length(dirGaze_resampled);
                    
                    % compute probabilities for each block
                    probFix_blocks.intervalNames = {'Preview','RT','Reach','Feedback'};
                    probFix_blocks.intervals = cumsum([0 nsNew]);
                    probFix_blocks.blockNames = ExpDetails(d).blockNames;
                    for b = 1 : nBlocks
                        probFix_blocks.(days{d}).(blocks{b}).start = NaN(sum(nsNew),1);
                        probFix_blocks.(days{d}).(blocks{b}).pearls = NaN(sum(nsNew),nPearlBins);
                        currBlock = b==blockNumber(:,d);
                        nAnalyzed_curr = sum(D_new.analyzedGazeData(currBlock));
                        dirGaze_curr = dirGaze_resampled(:,currBlock);
                        distGaze_curr = distGaze_resampled(:,currBlock);
                        probFix_blocks.(days{d}).(blocks{b}).start(1:nsTotal) = ...
                            sum(distGaze_curr<minMaxPercDistance(1),2)/nAnalyzed_curr;
                        for iPearl = 1 : nPearlBins
                            probFix_blocks.(days{d}).(blocks{b}).pearls(1:nsTotal,iPearl) = ...
                                sum(dirGaze_curr>pearlBinAngles(iPearl,1) & dirGaze_curr<pearlBinAngles(iPearl,2),2)/nAnalyzed_curr;
                        end
                    end
                    
                    % compute probabilities for each subblock in rotation block
                    probFix_subblocks.intervalNames = {'Preview','RT','Reach','Feedback'};
                    probFix_subblocks.intervals = cumsum([0 nsNew]);
                    for bsub = 1 : length(uniSubBlock_rot)
                        probFix_subblocks.(days{d}).(subBlocks_rot{bsub}).start = NaN(sum(nsNew),1);
                        probFix_subblocks.(days{d}).(subBlocks_rot{bsub}).pearls = NaN(sum(nsNew),nPearlBins);
                        currBlock = uniSubBlock_rot(bsub)==subBlockNumber;
                        %if bsub==1; currBlock(find(currBlock,8)) = false; end % exclude first 8 trials
                        nAnalyzed_curr = sum(D_new.analyzedGazeData(currBlock));
                        dirGaze_curr = dirGaze_resampled(:,currBlock);
                        distGaze_curr = distGaze_resampled(:,currBlock);
                        probFix_subblocks.(days{d}).(subBlocks_rot{bsub}).start(1:nsTotal) = ...
                            sum(distGaze_curr<minMaxPercDistance(1),2)/nAnalyzed_curr;
                        for iPearl = 1 : nPearlBins
                            probFix_subblocks.(days{d}).(subBlocks_rot{bsub}).pearls(1:nsTotal,iPearl) = ...
                                sum(dirGaze_curr>pearlBinAngles(iPearl,1) & dirGaze_curr<pearlBinAngles(iPearl,2),2)/nAnalyzed_curr;
                        end
                    end
                    
                end
                
            end % if ~all(noGazeData)
            
            %% Get timing and angles
            
            % preview, reaction, and movement time
            targetPreviewDur(1:n,d) = Exp.timing.previewDur;
            RT(1:n,d) = D_new.tTargetGoLeaveRingEnd(:,3)-D_new.tTargetGoLeaveRingEnd(:,2);
            MT(1:n,d) = D_new.tTargetGoLeaveRingEnd(:,4)-D_new.tTargetGoLeaveRingEnd(:,3);
            
            % hit angle, explicit angles and fixation angle
            cursorRotation(1:n,d) = Exp.cursorRotation;
            hitAngle_cursor(1:n,d) = D_new.hitAngle_cursor(:,1);
            hitAngle_hand(1:n,d) = D_new.hitAngle_hand(:,1);            
            explicitAngle(1:n,d) = D_new.explicit1Angle;
            explicitAngle_outliersRemoved(:,d) = explicitAngle(:,d);
            explicitAngle_outliersRemoved(explicitAngle(:,d)<-90 | explicitAngle(:,d)>45,d) = NaN; % remove outliers
            if ~all(noGazeData)
                fixAngle_preview_closestA(1:n,d) = fixAngles.(days{d}).preview.closestA;
                fixAngle_preview_closestOppA(1:n,d) = fixAngles.(days{d}).preview.closestOppA;
            else
                fixAngle_preview_closestA(1:n,d) = NaN;
                fixAngle_preview_closestOppA(1:n,d) = NaN;
            end
            
            fprintf('\n')
        end % end of loop over days
        
        %% Bin cursor and hand angle, explicit angle and fixation angle
        
        % remove report trials from RT and MT
        RT_all = RT;
        RT(~noReport) = NaN;
        MT_all = MT;
        MT(~noReport) = NaN;
        
        % angular error and report averaged in 8 trial bins
        nTrials = length(hitAngle_cursor);
        iBin = 1 : 8 : nTrials+1;
        for b = 1 : length(iBin)-1
            % cursor hit angle
            a = hitAngle_cursor(iBin(b):iBin(b+1)-1,:);
            mnBinnedHitAngle_cursor(b,:) = nanmean(a);
            medBinnedHitAngle_cursor(b,:) = nanmedian(a);
            sdBinnedHitAngle_cursor(b,:) = nanstd(a);
            % hand hit angle
            a = hitAngle_hand(iBin(b):iBin(b+1)-1,:);
            mnBinnedHitAngle_hand(b,:) = nanmean(a);
            medBinnedHitAngle_hand(b,:) = nanmedian(a);
            sdBinnedHitAngle_hand(b,:) = nanstd(a);
            % RT
            a = RT(iBin(b):iBin(b+1)-1,:);
            mnBinnedRT(b,:) = nanmean(a);
            medBinnedRT(b,:) = nanmedian(a);
            sdBinnedRT(b,:) = nanstd(a);
            % report angle
            a = explicitAngle(iBin(b):iBin(b+1)-1,:);
            mnBinnedExplicitAngle(b,:) = nanmean(a);
            medBinnedExplicitAngle(b,:) = nanmedian(a);
            sdBinnedExplicitAngle(b,:) = nanstd(a);
            a = explicitAngle_outliersRemoved(iBin(b):iBin(b+1)-1,:);
            mnBinnedExplicitAngle_outliersRemoved(b,:) = nanmean(a);
            sdBinnedExplicitAngle_outliersRemoved(b,:) = nanstd(a);
            % fixation angle closest to ideal aimpoint - during preview
            if any(blockNumber(iBin(b),:)>=rotBlockNumber) % rotation or washout block
                a = fixAngle_preview_closestA(iBin(b):iBin(b+1)-1,:);
                a(:,sum(~isnan(a))==1) = NaN;  % do not compute mean when there is only 1 value
                mnBinnedFixAngle_preview_closestA(b,:) = nanmean(a);
                medBinnedFixAngle_preview_closestA(b,:) = nanmedian(a);
                sdBinnedFixAngle_preview_closestA(b,:) = nanstd(a);
                a = fixAngle_preview_closestOppA(iBin(b):iBin(b+1)-1,:);
                a(:,sum(~isnan(a))==1) = NaN;  % do not compute mean when there is only 1 value
                mnBinnedFixAngle_preview_closestOppA(b,:) = nanmean(a);
                medBinnedFixAngle_preview_closestOppA(b,:) = nanmedian(a);
                sdBinnedFixAngle_preview_closestOppA(b,:) = nanstd(a);
            end
        end
        iBin = iBin(1:end-1)';
        
        %% Compute implicit angle
        
        % subtraction of explicit from hand angle
        mnBinnedImplicitAngle = mnBinnedHitAngle_hand - mnBinnedExplicitAngle;
        mnBinnedImplicitAngle_fix = mnBinnedHitAngle_hand - mnBinnedFixAngle_preview_closestA;
        if ~isempty(washoutBlockNumber)
            bin_WO = find(blockNumber(iBin,1)==washoutBlockNumber);
            mnBinnedImplicitAngle_fix(bin_WO,:) = mnBinnedHitAngle_hand(bin_WO,:) - mnBinnedFixAngle_preview_closestOppA(bin_WO,:);
        end
        
        %% Save individual results
        
        Results.blockNo = blockNumber;
        Results.trialNo = [trialNumber (1:length(trialNumber))'];
        Results.trialType = cell(nTrials,nDays);
        Results.trialType(~isnan(trialNumber) & noReport) = {'noReport'};
        Results.trialType(~isnan(trialNumber) & ~noReport) = {'report'};
        Results.cursorRotation = cursorRotation;
        Results.feedbackAndMTgood = feedbackAndMTGood;
        Results.targetPreviewDur = targetPreviewDur;
        Results.RT = RT;
        Results.RT_all = RT_all;
        Results.MT = MT;
        Results.MT_all = MT_all;
        Results.hitAngle_cursor = hitAngle_cursor;
        Results.hitAngle_hand = hitAngle_hand;
        Results.explicitAngle = explicitAngle;
        Results.explicitAngle_outliersRemoved = explicitAngle_outliersRemoved;
        Results.analyzedGazeData = analyzedGazeData;
        Results.fixAngle_closestA = fixAngle_preview_closestA;
        Results.fixAngle_closestOppA = fixAngle_preview_closestOppA;
        if any(analyzedGazeData)
            Results.fixAngles = fixAngles;
            Results.iFixTarget_aimFix_rot = iFixTarget_aimFix_rot;
            Results.percTimeFix = percTimeFix;
            Results.percTimeFix_blocks = percTimeFix_blocks;
            Results.probFix_blocks = probFix_blocks;
            Results.probFix_subblocks = probFix_subblocks;
        end
        % binned data
        Results.blockNo_bins = blockNumber(iBin,:);
        Results.trialNo_bins = [trialNumber(iBin,:) iBin];
        Results.mnBinnedHitAngle_cursor = mnBinnedHitAngle_cursor;
        Results.medBinnedHitAngle_cursor = medBinnedHitAngle_cursor;
        Results.sdBinnedHitAngle_cursor = sdBinnedHitAngle_cursor;
        Results.mnBinnedHitAngle_hand = mnBinnedHitAngle_hand;
        Results.medBinnedHitAngle_hand = medBinnedHitAngle_hand;
        Results.sdBinnedHitAngle_hand = sdBinnedHitAngle_hand;
        Results.mnBinnedRT = mnBinnedRT;
        Results.medBinnedRT = medBinnedRT;
        Results.sdBinnedRT = sdBinnedRT;
        Results.mnBinnedExplicitAngle = mnBinnedExplicitAngle;
        Results.medBinnedExplicitAngle = medBinnedExplicitAngle;
        Results.sdBinnedExplicitAngle = sdBinnedExplicitAngle;
        Results.mnBinnedExplicitAngle_outliersRemoved = mnBinnedExplicitAngle_outliersRemoved;
        Results.sdBinnedExplicitAngle_outliersRemoved = sdBinnedExplicitAngle_outliersRemoved;
        Results.mnBinnedImplicitAngle = mnBinnedImplicitAngle;
        Results.mnBinnedImplicitAngle_fix = mnBinnedImplicitAngle_fix;
        Results.mnBinnedFixAngle_closestA = mnBinnedFixAngle_preview_closestA;
        Results.medBinnedFixAngle_closestA = medBinnedFixAngle_preview_closestA;
        Results.sdBinnedFixAngle_closestA = sdBinnedFixAngle_preview_closestA;
        Results.mnBinnedFixAngle_closestOppA = mnBinnedFixAngle_preview_closestOppA;
        Results.medBinnedFixAngle_closestOppA = medBinnedFixAngle_preview_closestOppA;
        Results.sdBinnedFixAngle_closestOppA = sdBinnedFixAngle_preview_closestOppA;
        % criteria for reach and gaze analysis
        maxRT = unique(Exp.timing.maxRT);
        if length(maxRT)==1
            maxRT_str = num2str(Exp.timing.maxRT(1),'0-%.1f s');
        else
            maxRT_str = num2str(maxRT',['0-%.1f' repmat(' or 0-%.1f',1,length(maxRT)-1) ' s']);
        end
        if ~isfield(Exp.timing,'maxMT'); Exp.timing.maxMT = Exp.timing.maxReachDur; end % temp
        Results.criteria = {'For analysis of reach and gaze analysis: ';...
            ['FeedbackGood during testing: '];...
            ['RT ' maxRT_str ', max MT ' num2str(min([maxMT unique(Exp.timing.maxMT','stable')])) ' s'];...
            ['For analysis of fixations: '];...
            ['Min ' num2str(100-maxPercNoGaze) ' percent of gaze data (no blink/pupil lost) during preview and reach'];...
            ['Fixation distance between ' num2str(minMaxPercDistance) ' of reach distance'];...
            ['For analysis of reported values: '];...
            ['Classified as outlier if <-90 or >45 deg']};
        
        % check if file exists
        if nDays>1
            fileName_ = [Exp.subjName '.mat'];
            fileName = [Exp.subjFolder(1:end-5) '.mat'];
            if ~strcmp(fileName_,fileName) % temp
                disp('Check filename'); keyboard
            end
        end
        if exist([saveToPath fileName],'file') == 2 % check if file does not exist yet
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
        plotIndResultsVMRgaze(Results,ExpDetails,savePlots,saveFigsToPath)
    end
    
end % end of loop over subjects

% save copy of Matlab code used to process data
fileDate = datestr(now,'yyyymmmmdd');
filePath = mfilename('fullpath');
[fList,~] = matlab.codetools.requiredFilesAndProducts([filePath '.m']);
fList = [[filePath '.m'] fList];
fList = fList(cellfun(@isempty,strfind(fList,'ExpAlgemeen')));
for f = 1 : length(fList)
    mfilePath = fList{f};
    iSlash = strfind(mfilePath,'/');
    mfileName = mfilePath(iSlash(end)+1:end);
    mfileCopyName = sprintf('%s_copy_%s',mfileName,fileDate);
    copyfile(mfilePath,[saveToPath mfileCopyName]);
end

fprintf('\n'); disp('Done. Press continue to end.')
keyboard
close all
