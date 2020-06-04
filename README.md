# Analysis tools

## General
* **filterGazeData.m** - Remove blinks and filter raw gaze data from an Eyelink eye tracker.
* **findIntervals.m** - Extension of Matlab's ```find``` function. Find the start and end row indices of all intervals during which one or multiple conditions are nonzero for a minimum duration.
* **replaceErrorTrialsInStruct.m** - Replace values in a data struct with NaNs for trials categorized as error trials.
* **resampleGazeIntervals.m** - Resample gaze intervals to standard intervals for the purpose of averaging across trials. 
* **saccadeOnsetOffset2** - Find the onset and offset of saccades in eye movement data using a velocity threshold.
* **selectFiles.m** - Manually select files or folders matching a name (can include wildcards). 
* **xy2compassAngle.m** - Convert xy position to angle in degrees where 0 degrees is north and clockwise rotations are positive. 

## Tablet
Code for analyzing data collected on tablet setup
* **analyzeVMRgazeData.m** - Analyze gaze data of visuomotor rotation experiments
* **calcIndResultsVMRgaze.m** - Calculate individual participant results of visuomotor rotation experiments with gaze tracking.
* **createRLanimation.m** - Create avi animation of reward learning task
* **createVMRanimation.m** - Create avi animation of visuomotor rotation task
* **getColumnIndex.m** - Get index of data column by variable name
* **getHeaderValue.m** - Get value from text header of datafile by variable name
* **getTabletExpInfo.m** -  Get tablet experiment information from header of datafile
* **plotIndResultsVMRgaze.m** - Create plots with individual participant results calculated in ```calcIndResultsVMRgaze```
* **processRLtabletData.m** - Process raw data from reward learning experiments on tablet setup
* **processVMRtabletData.m** -  Process raw data from visuomotor rotation experiments on tablet setup
* **readFlanData.m** - Import Flanagan lab text data file
