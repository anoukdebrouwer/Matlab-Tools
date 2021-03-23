# Analysis tools

* **filterGazeData** - Remove blinks and filter raw gaze data from an Eyelink eye tracker.
* **findIntervals** - Extension of Matlab's ```find``` function. Find the start and end row indices of all intervals during which one or multiple conditions are nonzero for a minimum duration.
* **getFolderName** - Get the name of a folder at a certain level in the full path.
* **replaceErrorTrialsInStruct** - Replace values in a data struct with NaNs for trials categorized as error trials.
* **resampleGazeIntervals** - Resample gaze intervals to standard intervals for the purpose of averaging across trials. 
* **saccadeOnsetOffset2** - Find the onset and offset of saccades in eye movement data using a velocity threshold.
* **saveCopyOfCode** - Save copy of code (e.g., with processed data) for future reference.
* **selectFiles** - Manually select files or folders matching a name (can include wildcards). 
* **xy2compass** - Convert xy position to compass angle in degrees where 0 degrees is 'north' or up and clockwise rotations are positive. 
