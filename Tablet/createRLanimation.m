function createRLanimation(dataPath,trials)
% createRLanimation  Create avi animation of RL task

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

% define path where animation is saved
iSep = strfind(dataPath,filesep);
projectPath = dataPath(1:iSep(end-1)); % two levels up

%% Load data

load(dataPath)
t = trials(1);
b = Exp.blockNo(t);

%% Get stimuli
% to do: load raw data file and get sizes and colors from header

fs = 50;
monitorSize = Exp.setup.monitorSize_mm;
startXY = [0 0];
startRadius = Exp.stim.startRadius;
targetXY = Exp.targetXY(b,:);
targetSize = Exp.stim.targetWidthHeight;
cursorRadius = Exp.stim.cursorRadius;

% target path
targetPathCurvature = Exp.targetPathCurvature(b);
targetPathNumberOfCurves = Exp.targetPathNumberOfCurves(b);
targetPathVisible = Exp.targetPathVisible(b);
yTargetPath = (0:10:targetXY(2))';
y_norm = yTargetPath/targetXY(2);
x_norm = (0*y_norm) + (targetPathCurvature*sin(pi*targetPathNumberOfCurves*y_norm));
xTargetPath = x_norm*targetXY(2);

% visible path
visiblePathCurvature = Exp.visiblePathCurvature(b);
visiblePathNumberOfCurves = Exp.visiblePathNumberOfCurves(b);
visiblePathVisible = Exp.visiblePathVisible(b);
if isnan(visiblePathVisible)
    visiblePathVisible = 0;
end
yVisiblePath = (0:10:targetXY(2))';
y_norm = yVisiblePath/targetXY(2);
x_norm = (0*y_norm) + (visiblePathCurvature*sin(pi*visiblePathNumberOfCurves*y_norm));
xVisiblePath = x_norm*targetXY(2);

%% Open video file

v = VideoWriter([projectPath 'RLtrials_' Exp.subjName '.avi']);
v.FrameRate = fs;
open(v);
clear M
ni=0;

%% Create figure

figure; hold on
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) monitorSize*3]);
set(gca,'position',[0 0 1 1],'units','normalized')
rectangle('Position',[-monitorSize(1)/2 -monitorSize(2)/2+100 monitorSize],'FaceColor','k'); % background
axis equal; axis([-monitorSize(1)/2 monitorSize(1)/2 -monitorSize(2)/2+100 monitorSize(2)/2+100])
set(gca,'Xtick',[]); set(gca,'Ytick',[])

%% Get trial timing and cursor movement

for t = trials
    
    % get time and cursor movement
    time_tr_ = D.time{t};
    trialState_ = D.trialState{t};
    xyCursor_ = D.xyCursor{t};
    % define when cursor changed position
    [xyCursor_new,i_new] = unique(xyCursor_,'rows','stable');
    time_new = time_tr_(i_new)-time_tr_(i_new(1));
    trialState_new = trialState_(i_new);
    % resample
    time_r = (time_new(1):1/fs:time_new(end))';
    xyCursor_r = interp1(time_new,xyCursor_new,time_r);
    % get timing of events
    iHit = find(xyCursor_r(:,2)>=targetXY(2),1);
    
    % get score
    score = D.score(t);
    
    %% Create animation
    
    % plot border, start, path, target, and cursor
    brdr = rectangle('Position',[0-320/2 -25+100-170/2 320 170],'EdgeColor','w'); % border
    if targetPathVisible
        tpath = plot(xTargetPath,yTargetPath,'color',[0.5 0.5 0.5]); % rewarded path
    end
    if visiblePathVisible
        vpath = plot(xVisiblePath,yVisiblePath,'color',[0.5 0.5 0.5]); % visible path
    end
    strt = rectangle('Position',[0-startRadius 0-startRadius startRadius*2 startRadius*2],...        % start
        'Curvature',[1 1],'EdgeColor','none','FaceColor','w');
    trgt = rectangle('Position',[targetXY-targetSize/2 targetSize],'EdgeColor','r','FaceColor','r'); % target
    crsr = rectangle('Position',[0-cursorRadius 0-cursorRadius cursorRadius*2 cursorRadius*2],... % cursor
        'Curvature',[1 1],'EdgeColor','none','FaceColor',[0 1 1]);
    % start trial: reach
    i = 1;
    while i<iHit
        crsr.Position = [xyCursor_r(i,:)-startRadius startRadius*2 startRadius*2];
        M(i+ni) = getframe(gcf);
        writeVideo(v,M(i+ni))
        i = i+1;
    end
    % reach complete, delete cursor and show score
    crsr.delete
    txt = text(0,targetXY(2)-10,num2str(score),'HorizontalAlignment','center',...
        'FontSize',30,'FontWeight','bold','color','w');
    while i<(iHit+fs) % +1s
        M(i+ni) = getframe(gcf);
        writeVideo(v,M(i+ni))
        i = i+1;
    end
    % iti
    brdr.delete
    if exist('tpath','var'); tpath.delete; end
    if exist('vpath','var'); vpath.delete; end
    trgt.delete
    strt.delete
    crsr.delete
    txt.delete
    while i<(iHit+2*fs) % +1s
        M(i+ni) = getframe(gcf);
        writeVideo(v,M(i+ni))
        i = i+1;
    end
    M(i+ni) = getframe(gcf);
    writeVideo(v,M(i+ni))
    ni = ni+i;
end

close(v)
close(gcf)
