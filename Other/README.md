# Other Matlab tools

## findIntervals.m
Extension of Matlab's ```find``` function. Find the start and end row indices of all intervals during which one or multiple conditions (i.e., columns in the data  matrix) are nonzero. The minimum length of the interval (i.e., rows in the data matrix) during which the conditions need to be nonzero is specified by the user.

Example: 
```
%% find when both the x and y position of the measured hand position are in the target area for at least 10 measurement samples
xInTarget = x > (xtarget-0.5*targetWidth) & x < (xtarget+0.5*targetWidth);
yInTarget = y > (ytarget-0.5*targetHeight) & y < (ytarget+0.5*targetHeight);
[i,iStartEnd] = findIntervals([xInTarget yInTarget],10);
```

## selectFiles.m
Manually select files or folders matching a name (can include wildcards) from a list printed in the command window. Return a struct with file info.

Example:

<img src="/Images/selectFiles_example.png" width="600">

## xy2compassAngle.m
Convert xy position to angle in degrees where 0 degrees is north and clockwise rotations are positive. 

Example:
```
x = [0,10,0,-10];
y = [10,0,-10,0];
xy2compassAngle(x,y)
```
Returns ```0    90   180   270```

