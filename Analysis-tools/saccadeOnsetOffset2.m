function [onsets,offsets] = saccadeOnsetOffset2(xyvpGaze,searchWin,fs,sacDetParams)
% saccadeOnsetOffset2  Find the onset and offset of saccades using a
% velocity threshold
%
% [onsets,offsets] = saccadeOnsetOffset2(xyvpGaze,searchWin,fs) returns 
% a vector onsets with the onset indices and a vector offsets with the 
% offset indices of all saccades in the gaze data in xyvpGaze, within the 
% search window searchWin [iStart iEnd]. xyvpGaze should be a matrix that 
% contains the x and y gaze position and the resultant gaze velocity in 
% degrees of visual angle (filtered and blinks removed), and pupil area in 
% the columns. fs is the sampling frequency (e.g., 500) of these data.
% 
% The algorithm looks for a peaks in gaze velocity above velocity threshold 
% vTh for a minimum duration minPeakDur. It then searches backwards in
% time for the saccade onset, defined as the last sample below velocity
% threshold vOn for duration minOnDur, and forward in time for the saccade
% offset, defined as the first sample below velocity threshold vOn for
% duration minOnDur. Saccades with a displacement <minAmpl or >maxAmple 
% will be discarded. If parameter inclEdges is set to TRUE, velocity peaks 
% at the edges of the search window for which only onset or offset is found
% will be included. 
%
% If no saccade is detected, onsets and offsets will be 0 if there was no 
% peak in gaze velocity, or if a velocity peak was detected, but no onset 
% or offset were found, or gaze displacement was not sufficient. The 
% outputs will be -1 if there is no data during 95-100% of the search window.

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

plotOnOff = false; % set to true for visual inspection

% set default parameters
vTh = 75;           % velocity threshold for saccade detection
vOn = 30;           % velocity threshold for saccade onset and offset detection
vMax = 1000;        % maximum gaze velocity
minPeakDur = 4;     % minimum duration of velocity peak in ms (v>vTh)
minOnDur = 10;      % minimum duration of velocity<vOn for onset and offset in ms
minAmpl = 2;        % minimum saccade amplitude
maxAmpl = 100;      % maximum saccade amplitude
combineSacInt = 20; % maximum interval between offset and onset of subsequent
                    % saccades for which saccades will be combined
inclEdges = false;  % include velocity peaks at the edges of the search
                    % window as saccades, useful to filter out saccades
                    % for detecting fixations

% overwrite parameters if specified
if nargin == 4
    if isfield(sacDetParams,'vTh'); vTh = sacDetParams.vTh; end
    if isfield(sacDetParams,'vOn'); vOn = sacDetParams.vOn; end
    if isfield(sacDetParams,'vMax'); vMax = sacDetParams.vMax; end
    if isfield(sacDetParams,'minPeakDur'); minPeakDur = sacDetParams.minPeakDur; end
    if isfield(sacDetParams,'minOnDur'); minOnDur = sacDetParams.minOnDur; end
    if isfield(sacDetParams,'minAmpl'); minAmpl = sacDetParams.minAmpl; end
    if isfield(sacDetParams,'maxAmpl'); maxAmpl = sacDetParams.maxAmpl; end
    if isfield(sacDetParams,'combineSacInt'); combineSacInt = sacDetParams.combineSacInt; end
    if isfield(sacDetParams,'inclEdges'); inclEdges = sacDetParams.inclEdges; end
elseif nargin>4
    disp('Too many inputs. Use struct with saccade detection parameters.')
    keyboard
end

% select data
if searchWin(2)>length(xyvpGaze) % if data ends before search window, decrease size of window
    searchWin(2) = length(xyvpGaze);
end
searchWin = searchWin(1):searchWin(2);
x = xyvpGaze(searchWin,1);
y = xyvpGaze(searchWin,2);
v = abs(xyvpGaze(searchWin,3));
pa = xyvpGaze(searchWin,4);

T = 1000/fs; % sample time in ms
onOffInt = round(1000/T); % interval for onset and offset detection in samples

% preallocate
onsets = [];
offsets = [];

%%%% VELOCITY PEAKS %%%%
% find where v>vTh and v<maxVel within search window
vHigh = v>vTh & v<vMax;
peakStart = find(diff(vHigh)==1);
peakEnd = find(diff(vHigh)==-1);
if ~isempty(peakStart) || ~isempty(peakEnd) % velocity peak is found
    
    % include peaks at the edges of the search window
    % if inclEdges=false, saccades will be discarded later
    if isempty(peakStart) || isempty(peakEnd)
        peakStart(isempty(peakStart)) = find(pa,1);
        peakEnd(isempty(peakEnd)) = find(pa,1,'last');
    end
    if peakStart(1)>peakEnd(1)
        peakStart = [find(pa,1); peakStart];
    end
    if peakEnd(end)<peakStart(end)
        peakEnd = [peakEnd; find(pa,1,'last')];
    end
    velPeaks = [peakStart peakEnd];
    
    % duration of velocity peaks
    durPeaks = (velPeaks(:,2)-velPeaks(:,1))*T;
    durOK = durPeaks>=minPeakDur;
    velPeaks = velPeaks(durOK,:);
    nPeaks = size(velPeaks,1);
    
    % find onset and offset for each peak
    low = [];
    on = zeros(1,nPeaks);
    off = zeros(1,nPeaks);
    for i = 1:nPeaks
        
        %%%% SACCADE ONSET %%%%
        % look back in time to find saccade onset (last sample with v<vOn)
        st = max([velPeaks(i,1)-onOffInt 1]); % start of search interval (>0)
        v1 = v(st:velPeaks(i,1)+1);
        pa1 = pa(st:velPeaks(i,1)+1);
        vLow = v1>0 & v1<vOn & pa1>0;
        lowStart = find(diff(vLow)==1);
        lowEnd = find(diff(vLow)==-1);
        if vLow(1) && (isempty(lowStart) || lowStart(1)>lowEnd(1))
            lowStart = [1; lowStart]; % start of first low velocity interval is before search window
        end
        low = [lowStart lowEnd];
        if ~isempty(low)
            durLow = (low(:,end)-low(:,1))*T;
            low = low(durLow>=minOnDur,:);
            if ~isempty(low) % onset is found
                on(i) = low(end,2) + st -1;
            end
        end
        
        %%%% SACCADE OFFSET %%%%
        % look forward in time to find saccade offset (first sample with v<vOn)
        ende = min([velPeaks(i,2)+onOffInt length(v)]); % end of search interval (<length(v))
        v2 = v(velPeaks(i,2):ende);
        pa2 = pa(velPeaks(i,2):ende);
        vLow = v2>0 & v2<vOn & pa2>0;
        lowStart = find(diff(vLow)==1);
        lowEnd = find(diff(vLow)==-1);
        if vLow(end) && (isempty(lowEnd) || lowEnd(end)<lowStart(end))
            lowEnd = [lowEnd; length(v2)]; % end of low velocity interval is after search window
        end
        low = [lowStart lowEnd];
        if ~isempty(low)
            durLow = (low(:,2)-low(:,1))*T;
            low = low(durLow>=minOnDur,:);
            if ~isempty(low) % offset is found
                off(i) = low(1,1) + velPeaks(i,2);
            end
        end
    end
    
    % if velocity peak(s) is (are) found
    if nPeaks>0
        
        % if inclEdges==true, include saccades at the edges of the search window
        if inclEdges
            if on(1)==0 && off(1)>0
                on(1) = find(pa,1);
            end
            if off(end)==0 && on(end)>0
                off(end) = find(pa,1,'last');
            end
        end
        
        % remove saccades for which no onset or offset is found
        onAndOff = on>0 & off>0;
        on = on(onAndOff);
        off = off(onAndOff);
        
        % get unique onsets and offsets
        [uniOnOff,~] = unique([on(:) off(:)],'rows');
        on = uniOnOff(:,1);
        off = uniOnOff(:,2);
        
        % plot velocity and onsets and offsets
        if plotOnOff
            plot(v); axis([0 length(v) 0 vTh*10]);
            vertline(on,'g:'); vertline(off,'r:')
        end
        
        % combine saccades if <combineSacInt ms between subsequent saccades
        if combineSacInt>0 && length(on)>1
            int = (on(2:end)-off(1:end-1))*T;
            shortInt = find(int<combineSacInt,1);
            while ~isempty(shortInt)
                on(shortInt+1) = on(shortInt);
                on(shortInt) = [];
                off(shortInt) = [];
                if length(on)>1
                    int = (on(2:end)-off(1:end-1))*T;
                    shortInt = find(int<combineSacInt,1);
                else
                    shortInt = [];
                end
            end
        end
        
        % compute amplitude
        dx = x(off)-x(on);
        dy = y(off)-y(on);
        ampl = sqrt(dx.^2 + dy.^2);
        amplOK = ampl>minAmpl & ampl<maxAmpl;
        
        % plot onsets and offsets for valid saccades
        if plotOnOff
            hold on; horline([vTh vOn]);
            vertline(on(amplOK),'g-'); vertline(off(amplOK),'r-'); hold off;
            keyboard
        end
        
        % saccade onsets and offsets
        onsets = on(amplOK)+searchWin(1)-1; % saccade onset sample
        offsets = off(amplOK)+searchWin(1)-1; % saccade offset sample
    end
else
    % no velocity peak
    onsets = 0;
    offsets = 0;
end

%%%% CHECK WHY SACCADE WAS NOT DETECTED %%%%
if isempty(onsets)
    if plotOnOff
        plot(v); ylim([0 400]); horline(vTh);
    end
    if sum(pa==0)>(0.90*length(searchWin))
        % no data during >90% of search window, blink or pupil lost
        onsets = -1;
        offsets = -1;
        if plotOnOff
            title('No data during >90% of search window')
        end
    elseif all(~durOK)
        % no velocity peak with sufficient duration
        onsets = 0;
        offsets = 0;
        if plotOnOff
            title('No velocity peak with sufficient duration')
        end
    elseif isempty(low)
        % no onset or offset found
        onsets = 0;
        offsets = 0;
    elseif all(~amplOK)
        % no saccade with sufficient amplitude
        onsets = 0;
        offsets = 0;
        if plotOnOff
            title('No saccade with sufficient amplitude')
        end
    else
        onsets = 0;
        offsets = 0;
        keyboard
    end
    if plotOnOff
        waitforbuttonpress
    end
end