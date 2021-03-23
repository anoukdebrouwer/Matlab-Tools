## Analysis tools

* **filterGazeData** - Remove blinks and filter raw gaze data from an Eyelink eye tracker.
* **findIntervals** - Extension of Matlab's ```find``` function. Find the start and end row indices of all intervals during which one or multiple conditions are nonzero for a minimum duration.
* **getFolderName** - Get the name of a folder at a certain level in the full path.
* **replaceErrorTrialsInStruct** - Replace values in a data struct with NaNs for trials categorized as error trials.
* **resampleGazeIntervals** - Resample gaze intervals to standard intervals for the purpose of averaging across trials. 
* **saccadeOnsetOffset2** - Find the onset and offset of saccades in eye movement data using a velocity threshold.
* **saveCopyOfCode** - Save copy of code (e.g., with processed data) for future reference.
* **selectFiles** - Manually select files or folders matching a name (can include wildcards). 
* **xy2compass** - Convert xy position to compass angle in degrees where 0 degrees is 'north' or up and clockwise rotations are positive.

## Plotting tools

* **createColorGradient** - Create an RGB color gradient between two colors, with the start and end color and the number of colors in the gradient specified by the user.
* **errorb** - Draw error bars. By Jonathan C. Lansey, retrieved from [Matlab File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/27387-create-healthy-looking-error-bars).
* **getxPositionOfBars** - Get the x position of bars in a bar plot, e.g. for drawing error bars.
* **horline** - Draw a horizontal line at a specified location in the current figure.
* **plotCorrelationMatrix** - Visualize the correlation between each two variables in the data as a colored cell in a correlation matrix, using a red-white-blue color gradient. Display the correlation coefficient in each cell with asterisks indicating the significance.
    Example:
    ```
    load carbig
    X = [MPG Displacement Acceleration Horsepower Weight rand(length(MPG),1)];
    labels = {'MPG','Displacement','Acceleration','Horsepower','Weight','Random numbers'};

    figure
    plotCorrelationMatrix(X,labels);
    ```
    <img src="/Plotting-tools/Images/correlationMatrix_example.png" width="500">

* **plotMeansWithDataPoints** - Create a boxplot-like plot displaying the mean or median value as a line, the variance as a shaded area, and individual values as points. By default, the mean and standard error of the mean are plotted, but alternative measures for central tendency and variance can be specified by the user. Colors can be specified to group variables. 
    Example:
    ```
    load fisheriris
    setosa = ismember(species,'setosa');
    versicolor = ismember(species,'versicolor');
    virginica = ismember(species,'virginica');
    sepal_length = [meas(setosa,1) meas(versicolor,1) meas(virginica,1)];
    sepal_width = [meas(setosa,2) meas(versicolor,2) meas(virginica,2)];

    figure
    colors = [0.8 0 0; 0 0.8 0; 0 0 0.8];
    plotMeansWithDataPoints([sepal_length sepal_width],[],[],colors,{'setosa','versicolor','virginica'});
    set(gca,'XTick',[2,5])
    set(gca,'XTickLabel',{'Length','Width'})
    ```
    <img src="/Plotting-tools/Images/meansWithDataPoints_example.png" width="500">

* **rgb** - Get the rgb value of a color by name. By Kristjan Jonasson, retrieved form [Matlab File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/24497-rgb-triple-of-color-name-version-2).
* **saveFigAsPDF** - Save the current figure as a pdf file. The file name, path, and font size are specified by the user. 
* **scaledFigure** - Create a new figure window whose width and height are independently scaled by a factor relative to the default figure size.
* **scatterWithLinearFit** - Create a scatter plot with a least-squares fitted line and a text box displaying the number of valid data points, correlation coefficient and corresponding p-value.
    Example:
    ```
    load carbig
    scatterWithLinearFit(Horsepower,Acceleration)
    xlabel('Horsepower')
    ylabel('Acceleration')
    ```
    <img src="/Plotting-tools/Images/scatterWithLinearFit_example.png" width="500">

* **suplabel** - Place a title, xlabel, or ylabel on a group of subplots. By Ben Barrowes, retrieved from [Matlab File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/7772-suplabel).
* **vertline** - Draw a vertical line at a specified location in the current figure.
