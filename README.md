# Analysis tools

## General
* **findIntervals.m** - Extension of Matlab's ```find``` function. Find the start and end row indices of all intervals during which one or multiple conditions (i.e., columns in the data  matrix) are nonzero. The minimum length of the intervals (i.e., rows in the data matrix) can be specified by the user. Example use case: find intervals where both the measured x and y gaze position are in the target for at least 10 samples.
* **replaceErrorTrialsInStruct.m** - Replace values in a data struct with NaNs for trials categorized as error trials.
* **selectFiles.m** - Manually select files or folders matching a name (can include wildcards) from a list printed in the command window. 
* **xy2compassAngle.m** - Convert xy position to angle in degrees where 0 degrees is north and clockwise rotations are positive. 

## Tablet
* **createVMRanimation.m** - Create avi animation of visuomotor rotation task
* **getColumnIndex.m** - Get index of data column by variable name
* **getHeaderValue.m** - Get value from text header of datafile by variable name
* **getTabletExpInfo.m** -  Get tablet experiment information from header of datafile
* **processVMRtabletData.m** -  Process raw data from VMR experiments on tablet setup
* **readFlanData.m** - Import Flanagan lab text data file
