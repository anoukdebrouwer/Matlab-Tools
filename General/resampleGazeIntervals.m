function [gaze_resampled] = resampleGazeIntervals(gaze,oldIntervals,newIntervals)
% resampleGazeIntervals  Resamples gaze data to standard intervals

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

for i = 1 : length(newIntervals)
    if newIntervals(i)>0
        int_old = oldIntervals(i):oldIntervals(i+1)-1;
        if length(int_old)<=1; % create extra sample if interval contains a single sample
            int_old = [oldIntervals(i)-2 oldIntervals(i)];
        end
        int_new = linspace(int_old(1),int_old(end),newIntervals(i));
        % resample
        gaze_resampled{i}(1:newIntervals(i)) = interp1(int_old,gaze(int_old),int_new);        
    end
end
