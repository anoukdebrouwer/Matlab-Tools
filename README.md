# Matlab-tools

## createColorGradient.m
Create an RGB color gradient between two colors, with the start and end color and the number of colors in the gradient specified by the user.

Example:
```
createColorGradient([0 1 0],[0 0 1],20)
```

## findIntervals.m
Find the start and end indices of all intervals during which one or multiple conditions (i.e., columns in the data  matrix) are nonzero. The minimum length of the interval (i.e., rows in the data matrix) during which the conditions need to be nonzero is specified by the user.

Example: 
```
%% find when both the x and y position of the measured hand position are in the target area for at least 10 measurement samples
xInTarget = x > (xtarget-0.5*targetWidth) & x < (xtarget+0.5*targetWidth);
yInTarget = y > (ytarget-0.5*targetHeight) & y < (ytarget+0.5*targetHeight);
[i,iStartEnd] = findIntervals([xInTarget yInTarget],10);
```

## plotCorrelationMatrix.m
Visualize the correlation between each two variables in the data as a colored cell in a correlation matrix, using a red-white-blue color gradient. Display the correlation coefficient in each cell with asterisks indicating the significance.

Example:
```
load carbig
X = [MPG Displacement Acceleration Horsepower Weight rand(length(MPG),1)];
labels = {'MPG','Displacement','Acceleration','Horsepower','Weight','Random numbers'};

figure
plotCorrelationMatrix(X,labels);
```
![](/Images/correlationMatrix_example.png)

## plotMeansWithDataPoints.m
Create a boxplot-like plot displaying the mean or median value as a line, the variance as a shaded area, and individual values as points. By default, the mean and standard error of the mean are plotted, but alternative measures for central tendency and variance can be specified by the user. Colors can be specified to group variables. 

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
![](/Images/meansWithDataPoints_example.png)

## saveFigAsPDF.m
Save the current figure as a pdf file. The file name, path, and font size are specified by the user. 

Example:
```
saveFigAsPDF('meansWithDataPoints_example',10)
```

## scaledFigure.m
Create a new figure window whose width and height are independently scaled by a factor relative to the default figure size.

Example:
```
scaledFigure(2,1)
```

## scatterWithLinearFit.m
Create a scatter plot with a least-squares fitted line and a text box displaying the number of valid data points, correlation coefficient and corresponding p-value

Example:
```
load carbig
scatterWithLinearFit(Horsepower,Acceleration)
xlabel('Horsepower')
ylabel('Acceleration')
```
![](/Images/scatterWithLinearFit_example.png)

## selectFiles.m
Manually select files or folders matching a name (can include wildcards) from a list printed in the command window. Return a struct with file info.

Example:
![](/Images/selectFiles_example.png)
