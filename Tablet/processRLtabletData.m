function processRLtabletData(projectPath,plotTrials)
% processRLtabletData  Process raw data from RL experiments on tablet setup.
%
% processRLtabletData loads, processes, and visualizes hand movement data
% collected in a reward-based motor learning (RL) experiment on a tablet
% setup at Queen's University:
% - Wacom tablet with pen in Abramsky Hall
% - Wacom tablet in mock MRI scanner in Craine Building
% - custom MRI-compatible tablet in fMRI scanner
%
% processRLtabletData(projectPath) asks the user to specify a project
% folder (if unspecified) and select one or multiple participants or
% sessions. The raw data should be stored in a subfolder named '1_RawData',
% with a separate folder containing the raw data (one file per trial) of
% each participant and session. The function loads and cleans the data of
% each trial and saves a single file per participant and session,
% containing a struct 'Exp' with information about the experiment and
% trials, and a struct 'D' with data. The processed data will be saved in a
% subfolder named '2_ProcessedData'.
%
% processRLTabletData(projectPath,TRUE) plots each trial separately for
% visual inspection during processing, and waits for the user to click a
% button to continue to the next trial. The plot is an x y plot of the
% target path that the participant had to copy or find, the participant's
% cursor path, and the score shown to the participant.

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

close all;

if nargin==0
    projectPath = [];
    plotTrials = false;
end

% set project data path
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
dataPath = [projectPath '/1_RawData/'];
saveToPath = [projectPath '/2_ProcessedData/'];
if ~exist(saveToPath,'dir')
    mkdir(saveToPath)
end

% select subjects
cd(dataPath)
subjFolders = selectFiles('*','folders');
nSubj = length(subjFolders);

% open figure
if plotTrials
    fig1 = scaledFigure(1.5,1);
    colors = get(gca,'colororder');
    fadedcolors = colors+(1-colors)*0.7;
end

%% Loop over subjects

for s = 1 : nSubj
    
    % create structs for saving
    Exp=[];
    D=[];
    
    % get data files
    fprintf('\nLoading %s ...',subjFolders(s).name)
    subjDir = [dataPath '/' subjFolders(s).name '/'];
    controlFiles = dir([subjDir '*.txt']);
    Exp.controlFiles = {controlFiles.name}';
    fprintf(['\nControl file(s): ' repmat('%s,',1,length(controlFiles)-1) '%s\n'],...
        Exp.controlFiles{:})
    dataFiles_unsorted = dir([subjDir '*RC*.dat']);
    % sort trials in order of presentation,
    % and remove trials that were replayed because of incorrect timing
    [~,order] = sort([dataFiles_unsorted.datenum]);
    dataFiles_sorted = dataFiles_unsorted(order,:);
    replayedTrials = cellfun(@(x) ~isempty(regexp(x,'\w*_R0\d.dat','once')),...
        {dataFiles_sorted.name}');
    dataFiles = dataFiles_sorted(~replayedTrials,:);
    fprintf('\n%d trials replayed because of incorrect timing\n',sum(replayedTrials))
    nTrials = length(dataFiles);
    [i1,i2] = regexp(dataFiles(end).name,regexptranslate('wildcard','B*T'));
    nBlocks = str2double(dataFiles(end).name(i1+1:i2-1));
    if nTrials==0
        disp('No data files found')
        keyboard
    end
    
    % preallocate variables - trials
    dateTime        = NaN(nTrials,6);
    trialNo         = NaN(nTrials,1);
    blockNo         = NaN(nTrials,1);
    trialType       = cell(nTrials,1);
    feedbackMessage = cell(nTrials,1);
    iTargetGoLeaveCompleteEnd = NaN(nTrials,5);
    tTargetGoLeaveCompleteEnd = NaN(nTrials,5);
    score           = NaN(nTrials,1);
    pathTotalErr    = NaN(nTrials,1);
    score_offline   = NaN(nTrials,1);
    pathTotalDiff   = NaN(nTrials,1);
    pathChange      = NaN(nTrials,1);
    
    % preallocate variables - blocks
    visiblePathCurvature        = NaN(nBlocks,1);
    visiblePathNumberOfCurves   = NaN(nBlocks,1);
    visiblePathVisible          = NaN(nBlocks,1);
    targetPathCurvature         = NaN(nBlocks,1);
    targetPathNumberOfCurves    = NaN(nBlocks,1);
    targetPathVisible           = NaN(nBlocks,1);
    showScore       = NaN(nBlocks,1);
    highScore       = NaN(nBlocks,1);
    
    %% Loop over trials
    
    for t = 1 : nTrials
        
        % read in data to T structure
        T = ReadFlanData([subjDir dataFiles(t).name]);
        textdata = T.textdata(:,1);
        
        %% Get general experiment information
        
        % get general info from datafile of first trial
        if t==1
            
            Exp.expName = expName;
            Exp = getTabletExpInfo(Exp,textdata);
            Exp.subjFolder = subjFolders(s).name;
            Exp.rawdataFiles = rmfield(dataFiles,{'bytes','isdir','datenum'});
            
            % get number of target object
            objectText = [];
            lineNo = 0;
            while isempty(objectText)
                lineNo = lineNo+1;
                objectText = regexp(textdata{lineNo},'Object\w*Name\sTarget','match');
            end
            i = find(isstrprop(objectText{:},'digit'),1); % object number
            targetObjectNo_str = objectText{:}(i);
            
            % get start position (same for all trials)
            startXY_raw = [Exp.stim.startXY(1) Exp.stim.startXY(2)*(-1)];
            Exp.stim.startXY_raw = startXY_raw;
            Exp.stim.startXY = [0 0]; % data will be aligned to start
            
        end
        
        %% Get trial information
        
        % get trial information
        dateTimeStr = getHeaderValue(textdata,'Created');
        dateTime(t,:) = datevec(dateTimeStr(2:end-1),'ddd mmm dd HH:MM:SS yyyy');
        blockNo(t) = getHeaderValue(textdata,'BlockNumber');
        if t==1 || blockNo(t)>blockNo(t-1)
            firstTrialInBlock = true;
            trialNo(t) = 1; % trial number, restart numbering in new block
            b = blockNo(t);
        else
            firstTrialInBlock = false;
            trialNo(t) = trialNo(t-1)+1;
        end
        trialType{t} = getHeaderValue(textdata,'TrialType');
        feedbackTrial = strcmp(trialType{t},'ReachPathB'); % feedback trial
        abortedTrial = getHeaderValue(textdata,'TrialAborted'); % aborted trial
        if abortedTrial
            keyboard
            feedbackMessage{t} = 'Aborted';
        end
        
        % if trial was feedback trial (showing high score at end of block),
        % continue to next trial
        if feedbackTrial
            highScore(b) = getHeaderValue(textdata,'HighScore');
            continue
        end
        
        %% Get stimulus information per block
        
        if firstTrialInBlock
            
            % get target position
            targetXY_raw(1) = getHeaderValue(textdata,['Object' targetObjectNo_str 'PosX']);
            targetXY_raw(2) = getHeaderValue(textdata,['Object' targetObjectNo_str 'PosY'])*(-1); % default: down is positive
            Exp.stim.targetXY = targetXY_raw-startXY_raw; % target position relative to start position
            targetDistance = Exp.stim.targetDistance;
            
            % get target path specs - for score
            targetPathCurvature(b) = getHeaderValue(textdata,'CursorPathCurvature');
            targetPathNumberOfCurves(b) = getHeaderValue(textdata,'CursorPathGain');
            targetPathVisible(b) = getHeaderValue(textdata,'CursorPathVisible');
            
            % get visible path specs - for tracing or baseline
            visiblePathCurvature(b) = getHeaderValue(textdata,'VisibleCursorPathCurvature');
            visiblePathNumberOfCurves(b) = getHeaderValue(textdata,'VisibleCursorPathGain');
            visiblePathVisible(b) = getHeaderValue(textdata,'VisibleCursorPathVisible');
            
            % score at end of trial
            showScore(b) = getHeaderValue(textdata,'ShowFeedbackScore');
            
            %% Get path coordinates
            
            % get coordinates of target path for calculation of score
            yTargetPath = (0:10:targetDistance)';
            y_norm = yTargetPath/targetDistance;
            x_norm = targetPathCurvature(b)*sin(pi*targetPathNumberOfCurves(b)*y_norm);
            xTargetPath = x_norm*targetDistance;
            xTargetPath_score(:,b) = xTargetPath;
            yTargetPath_score(:,b) = yTargetPath;
            
            % get coordinates of path for scaling of score
            scaling = getHeaderValue(textdata,'CursorPathMaxCurvature');
            xMax_norm = scaling*sin(pi*targetPathNumberOfCurves(b)*y_norm);
            xMaxTargetPath = xMax_norm*targetDistance;
            
            % get coordinates of path for scaling of path changes
            xAmplOne_norm = 1*sin(pi*targetPathNumberOfCurves(b)*y_norm);
            xAmplOneTargetPath = xAmplOne_norm*targetDistance;
            
            % get coordinates of target path for plotting (finely sampled)
            yTargetPath = (0:1:targetDistance)';
            y_norm = yTargetPath/targetDistance;
            x_norm = targetPathCurvature(b)*sin(pi*targetPathNumberOfCurves(b)*y_norm);
            xTargetPath = x_norm*targetDistance;
            xTargetPath_plot(:,b) = xTargetPath;
            yTargetPath_plot(:,b) = yTargetPath;
            xyTargetPath_plot = [xTargetPath yTargetPath];
            
            % get coordinates of visible path for plotting (finely sampled)
            if targetPathVisible(b)
                xyVisiblePath_plot = xyTargetPath_plot;
            else
                yVisiblePath = (0:1:targetDistance)';
                y_norm = yVisiblePath/targetDistance;
                x_norm = visiblePathCurvature(b)*sin(pi*visiblePathNumberOfCurves(b)*y_norm);
                xVisiblePath = x_norm*targetDistance;
                xVisiblePath_plot(:,b) = xVisiblePath;
                yVisiblePath_plot(:,b) = yVisiblePath;
                xyVisiblePath_plot = [xVisiblePath yVisiblePath];
                if isnan(visiblePathCurvature(b))
                    xyVisiblePath_plot = xyVisiblePath_plot*NaN;
                end
            end
            
        end
        
        %% Get trial feedback
        
        % feedback message
        feedbackMessage{t} = getHeaderValue(textdata,'FeedbackState');
        
        % check if timeout occurred or feedbackmessage is missing
        col = getColumnIndex(T.colheaders,'TrialState');
        trialState = T.data(:,col);
        if any(trialState==14) || all(isnan(feedbackMessage{t}))
            % feedbackmessage is missing in datafile when timeout occured during trial state 1 or 2
            if isnan(feedbackMessage{t})
                str = sprintf('No feedback message in block %d trial %d',blockNo(t),trialNo(t));
                disp(str);
            end
            if any(trialState==13) % reach complete before timeout
                feedbackMessage{t} = 'Good';
            elseif ~any(trialState==13) && T.data(end,1)>(maxTrialDuration-0.002) % reach not complete before timeout
                feedbackMessage{t} = 'TimeOut';
            end
        end
        % if missing, add feedback message manually
        while isnan(feedbackMessage{t})
            disp('Add feedback message manually - check T.data'); keyboard
            % feedbackMessage{t} = 'Good';
        end
        
        % if trial timed out, move on to next trial
        if strcmp(feedbackMessage{t},'TimeOut')
            continue
        end
        
        %% Get raw data: time, trial states, pen, cursor, and score
        
        % remove data recorded before any target is visible (align to display delay)
        iStart = find(trialState==2,1); % start of recording
        if isempty(iStart); iStart=1; end
        data = T.data(iStart:end,:);
        
        % remove duplicate samples
        col = getColumnIndex(T.colheaders,'Time');
        time_tr = data(:,col);
        notNew = [0; diff(time_tr)==0];
        data = data(~notNew,:);
        
        % preallocate/clear variables
        nSamples = size(data,1);
        xyPen_raw = NaN(nSamples,2);
        xyCursor_raw = NaN(nSamples,2);
        xyStart_raw = ones(nSamples,1)*startXY_raw;
        
        % get time
        col = getColumnIndex(T.colheaders,'Time');
        time_tr = data(:,col);
        
        % get trial state names and durations
        col = getColumnIndex(T.colheaders,'TrialState');
        trialState = data(:,col);
        trialStateNames = Exp.timing.trialStateNames;
        trialStateDur = NaN(length(trialStateNames),1);
        for i = 1 : length(trialStateNames)
            iStateOn = find(trialState==i,1);
            iStateOff = find(trialState==i,1,'last');
            if ~isempty(iStateOn)
                trialStateDur(i) = time_tr(iStateOff) - time_tr(iStateOn) + 1/Exp.setup.fs;
            end
        end
        
        % get timing of events (from trial states)
        iTarget = find(trialState>=2,1);        % target is shown at end of holdstart
        iGo = iTarget;                          % go cue is simultaneous with target appearance
        iLeave = find(trialState==12,1);        % detection of movement start (online, when cursor reaches LeaveStartDistance)
        iComplete = find(trialState==13,1);     % reach completed (online, when cursor reaches target distance)
        iEnd = length(trialState);              % end of recording
        
        % get pen on tablet position
        col = getColumnIndex(T.colheaders,'TabletX');
        xyPen_raw(:,1) = data(:,col);
        col = getColumnIndex(T.colheaders,'TabletY');
        xyPen_raw(:,2) = data(:,col)*(-1);
        xyPen = xyPen_raw - xyStart_raw; % relative to start position
        
        % get cursor on screen position
        col = getColumnIndex(T.colheaders,'CursorPosX');
        xyCursor_raw(:,1) = data(:,col);
        col = getColumnIndex(T.colheaders,'CursorPosY');
        xyCursor_raw(:,2) = data(:,col)*(-1);
        xyCursor = xyCursor_raw - xyStart_raw; % relative to start position
        dCursor = sqrt(xyCursor(:,1).^2+xyCursor(:,2).^2); % distance from start
        
        % get start and end of reach from cursor position
        iStart = find(dCursor>Exp.timing.leaveStartDistance,1); % timing of leave start
        iHit = find(xyCursor(:,2)>targetDistance,1); % timing of hit
        
        % get score presented to participant
        score(t) = getHeaderValue(textdata,'Score');
        if showScore(b)
            score_shown = score(t);
        else
            score_shown = NaN;
        end
        
        %% Check for sampling issues (MRI-compatible tablet)
        % on MRI-compatible tablet, sampling issues can occur when the
        % finger doesn't sufficiently press down on the tablet,
        % classify as sampling issue if no samples were collected during reach
        
        % get number of samples during the reach
        xyCursor_reach = xyCursor(trialState==12,:);
        nSamples_reach = size(unique(xyCursor_reach,'rows'),1);
        % get the last sample before the cursor reaches the target distance
        xyCursor_beforeHit = xyCursor(1:iHit-1,:);
        nSamples_beforeHit = size(unique(xyCursor_beforeHit,'rows'),1);
        if nSamples_beforeHit>0
            iHit0 = find(xyCursor(:,2)==xyCursor_beforeHit(end,2),1); % last new sample before hit
        end
        % check number of samples and last cursor position before hit
        if nSamples_reach==0 || dCursor(iHit0)<Exp.timing.leaveStartDistance
            feedbackMessage{t} = 'SamplingIssue';
            if plotTrials
                plotCursorTrajectory(fig1,Exp,xyTargetPath_plot,xyVisiblePath_plot,...
                    xyCursor,[NaN NaN],score_shown)
                str = sprintf('Sampling issue: %s - block %d trial %d',...
                    subjFolders(s).name, blockNo(t), trialNo(t));
                title(str);
                keyboard
            end
        end
        
        % if trial is invalid, move on to next trial
        feedbackGood = strcmp(feedbackMessage{t},'Good');
        if ~feedbackGood
            score(t) = NaN;
            continue
        end
        
        % save timing
        iTargetGoLeaveCompleteEnd(t,:) = [iTarget iGo iStart iHit iEnd];
        tTargetGoLeaveCompleteEnd(t,:) = time_tr(iTargetGoLeaveCompleteEnd(t,:));
        
        %% Check score and calculate trial-by-trial path changes
        % the score is calculated by computing the distance between the
        % cursor position and the target path at every cm in the y direction
        
        % target path
        xTargetPath = xTargetPath_score(:,b);
        yTargetPath = yTargetPath_score(:,b);
        
        % store previous trial cursor positions to calculate path change
        xCursor_score_prev = NaN(length(yTargetPath),1);
        if ~firstTrialInBlock
            xCursor_score_prev = xCursor_score;
        end
        
        % interpolate cursor trajectory at every cm
        xCursor_score = zeros(length(yTargetPath),1);
        for i = 2 : length(yTargetPath)
            iy2 = find(xyCursor(:,2)>yTargetPath(i),1);
            iy1 = iy2-1;
            yRatio = (yTargetPath(i)-xyCursor(iy1,2)) / (xyCursor(iy2,2)-xyCursor(iy1,2));
            dx1 = yRatio * (xyCursor(iy2,1)-xyCursor(iy1,1));
            xCursor_score(i) = xyCursor(iy1,1) + dx1;
        end
        
        % calculate error and score
        maxTotalErr = sum(abs(xTargetPath)) + sum(abs(xMaxTargetPath));
        pathErr = xCursor_score-xTargetPath;
        pathTotalErr(t) = sum(abs(pathErr));
        score_offline(t) = round((maxTotalErr-pathTotalErr(t))/maxTotalErr*100);
        if score_offline(t)<0 % don't allow negative scores
            score_offline(t) = 0;
        end
        
        % compute change relative to previous trial
        if ~firstTrialInBlock
            pathDiff = xCursor_score-xCursor_score_prev;
            pathTotalDiff(t) = sum(abs(pathDiff));
            maxDiff = sum(abs(xAmplOneTargetPath)); % diff between straight path and ampl 1 curve
            pathChange(t) = 100-(round((maxDiff-pathTotalDiff(t))/maxDiff*100));
        end
        
        %% Store time, trial state, pen and cursor position in struct
        
        D.time{t,1}         = time_tr;
        D.trialState{t,1}   = trialState;
        D.xyPen{t,1}        = xyPen;
        D.xyCursor{t,1}     = xyCursor;
        
        %% Plot trials for visual inspection
        
        if plotTrials
            xyCursor_score = [xCursor_score,yTargetPath];
            plotCursorTrajectory(fig1,Exp,xyTargetPath_plot,xyVisiblePath_plot,...
                xyCursor,xyCursor_score,score_shown)
            hold on
            % add previous trial and path change
            % hold on
            % plot(xCursor_score_prev,yTargetPath,'.-','color',fadedcolors(1,:))
            % str = ['\Delta=' num2str(pathChange(t))];
            % text(xCursor_score_prev(end-1),yTargetPath(end-2),str,...
            %     'Color',fadedcolors(1,:),'HorizontalAlignment','center');
            % wait for button press
            title(['Block ' num2str(blockNo(t)) ' trial ' num2str(trialNo(t))])
            waitforbuttonpress % click/press key to continue
        end
        
    end % end of loop over trials
    
    % save scores of good trials and all trials
    score = repmat(score,1,2);
    score(strcmp(feedbackMessage,'SamplingIssue'),1) = NaN;
    
    %% Print number of good trials and error trials in command window
    
    feedbackTrial = strcmp(trialType,'ReachPathB');
    feedbackMessage_ = feedbackMessage(~feedbackTrial);
    
    nGood = sum(strcmp(feedbackMessage_,'Good'));
    nTimeOut = sum(strcmp(feedbackMessage_,'TimeOut'));
    nSamplingIssues = sum(strcmp(feedbackMessage_,'SamplingIssue'));
    str = sprintf('%d / %d trials with good timing',nGood,length(feedbackMessage_)); disp(str)
    str = sprintf('%d trials with timeout',nTimeOut); disp(str)
    str = sprintf('%d trials with sampling issues',nSamplingIssues); disp(str)
    
    %% Create structs for saving
    
    % create experiment struct
    Exp.processingDate = date;
    Exp.processingCode = 'processRLtabletData.m by Anouk de Brouwer';
    
    % stimulus information
    % trials - remove feedback trials
    Exp.blockNo = blockNo(~feedbackTrial);
    Exp.trialNo = trialNo(~feedbackTrial);
    Exp.feedbackBlockTrialNo = [blockNo(feedbackTrial) trialNo(feedbackTrial)];
    Exp.trialType = trialType(~feedbackTrial);
    % blocks
    Exp.scoreFeedback = showScore;
    Exp.targetXY = repmat(Exp.stim.targetXY,nBlocks,1);
    Exp.targetDistance = repmat(targetDistance,nBlocks,1);
    Exp.targetPathCurvature = targetPathCurvature;
    Exp.targetPathNumberOfCurves = targetPathNumberOfCurves;
    Exp.targetPathVisible = targetPathVisible;
    Exp.xTargetPath_plot = xTargetPath_plot;
    Exp.yTargetPath_plot = yTargetPath_plot;
    Exp.visiblePathCurvature = visiblePathCurvature;
    Exp.visiblePathNumberOfCurves = visiblePathNumberOfCurves;
    Exp.visiblePathVisible = visiblePathVisible;
    Exp.xVisiblePath_plot = xVisiblePath_plot;
    Exp.yVisiblePath_plot = yVisiblePath_plot;
    
    % data: feedback and timing
    D.time          = D.time(~feedbackTrial);
    D.trialState    = D.trialState(~feedbackTrial);
    D.feedbackMessage = feedbackMessage(~feedbackTrial);
    D.iTargetGoLeaveCompleteEnd = iTargetGoLeaveCompleteEnd(~feedbackTrial,:);
    D.tTargetGoLeaveCompleteEnd = tTargetGoLeaveCompleteEnd(~feedbackTrial,:);
    
    % data: hand movement
    % defined in trial loop:
    D.xyPen     = D.xyPen(~feedbackTrial);
    D.xyCursor  = D.xyCursor(~feedbackTrial);
    D.score     = score(~feedbackTrial);
    D.pathChange = pathChange(~feedbackTrial);
    D.description = {'All data is aligned to start of recording,';...
        'occurring after the hand is at start and before the stimuli are presented.';...
        'XY positions are in mm, time is in s.'};
    
    %% Save
    
    % check if file exists
    fileName = [subjFolders(s).name '.mat'];
    if exist([saveToPath fileName],'file') == 2 % check if file does not exist yet
        disp(['A file named ' fileName ' already exists.'])
        overwrite = input('Do you want to overwrite it? Yes(1) or no(0): ');
    else
        overwrite = 1;
    end
    
    % save file
    if overwrite == 1
        save([saveToPath fileName],'Exp','D');
        disp(['Saved ' saveToPath fileName])
    else
        disp('Data has not been saved')
    end
    
    if plotTrials==1
        disp('Click continue to go to next subject')
        keyboard
    end
    
end % end of loop over subjects

%% Save copy of code used to process data

mFilePath = mfilename('fullpath');
saveCopyOfCode(mFilePath,saveToPath)

end

function plotCursorTrajectory(fig,Exp,xyTargetPath,xyVisiblePath,xyCursor,xyCursor_score,score)

figure(fig); clf
colors = get(gca,'colororder');
% plot start and target
d = Exp.timing.leaveStartDistance;
targetXY = Exp.stim.targetXY;
plot(0,0,'ko'); hold on; plot(d*cosd(0:180),d*sind(0:180),'k:');
rectangle('Position',[targetXY-0.5*Exp.stim.targetWidthHeight Exp.stim.targetWidthHeight]);
% plot invisible and visible target path
ip = plot(xyTargetPath(:,1),xyTargetPath(:,2),'k--','linewidth',2);
vp = plot(xyVisiblePath(:,1),xyVisiblePath(:,2),'k-','linewidth',2);
% plot cursor
c = plot(xyCursor(:,1),xyCursor(:,2),'.-','color',colors(1,:));
% plot cursor locations for score
cs = plot(xyCursor_score(:,1),xyCursor_score(:,2),'x','color',c.Color);
% add score
text(targetXY(1),targetXY(2)-5,num2str(score),'HorizontalAlignment','Center');
% add legend and labels
targetDist = Exp.stim.targetDistance;
axis([-1.05*targetDist 1.05*targetDist -0.05*targetDist 1.05*targetDist])
daspect([1 1 1]);
if isnan(xyVisiblePath)
    legend([ip c cs],{'target path','cursor path','positions for score'},...
        'location','northwest');
else
    legend([ip vp c cs],{'target path','visible path','cursor path','positions for score'},...
        'location','northwest');
end
legend('boxoff')
xlabel('X (mm)'); ylabel('Y (mm)')

end