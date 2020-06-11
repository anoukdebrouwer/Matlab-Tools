# Plotting tools

## createColorGradient.m
Create an RGB color gradient between two colors, with the start and end color and the number of colors in the gradient specified by the user.

## errorb.m
Draw error bars. By Jonathan C. Lansey, retrieved from [Matlab File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/27387-create-healthy-looking-error-bars).

## getxPositionOfBars.m
Get the x position of bars in a bar plot, e.g. for drawing error bars.

## horline.m
Draw a horizontal line in the current figure.

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
<img src="/Images/correlationMatrix_example.png" width="500">

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
<img src="/Images/meansWithDataPoints_example.png" width="500">

## saveFigAsPDF.m
Save the current figure as a pdf file. The file name, path, and font size are specified by the user. 

## scaledFigure.m
Create a new figure window whose width and height are independently scaled by a factor relative to the default figure size.

## scatterWithLinearFit.m
Create a scatter plot with a least-squares fitted line and a text box displaying the number of valid data points, correlation coefficient and corresponding p-value

Example:
```
load carbig
scatterWithLinearFit(Horsepower,Acceleration)
xlabel('Horsepower')
ylabel('Acceleration')
```
<img src="/Images/scatterWithLinearFit_example.png" width="500">

## suplabel.m
Place a title, xlabel, or ylabel on a group of subplots. By Ben Barrowes, retrieved from [Matlab File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/7772-suplabel).

## vertline.m
Draw a vertical line in the current figure.
