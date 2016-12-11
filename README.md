# Categorical Scatter plot (Box Whisker plots with data points)
A replacement for the traditional box and whisker plots provided in MATLAB (command boxplot). The categorical scatter plots additionally shows the data points, which is useful the visualize the underlying distribution (similar to violin plots).

## Examples
The code is designed to be an extremely customizable alternate for the built in boxplot function in MATLAB. The syntax is very similar to that of boxplot. Below are some example usages for the function.

### Example 1: Using CategoricalScatterplot instead of boxplots
```
% load a sample dataset
load carsmall
% Use matlab's boxplot
boxplot(MPG,Origin)
title('Miles per Gallon by Vehicle Origin')
xlabel('Country of Origin')
ylabel('Miles per Gallon (MPG)')
```
The traditional box plot of the car sample dataset looks like this:
![](https://github.com/AbstractGeek/CategoricalScatterplot/img/box-plot.png "box plot")

```
% Use CategoricalScatterplot
CategoricalScatterplot(MPG,Origin)
title('Miles per Gallon by Vehicle Origin')
xlabel('Country of Origin')
ylabel('Miles per Gallon (MPG)')
```
The categoricalscatterplot of the same sample dataset looks like this:
![](https://github.com/AbstractGeek/CategoricalScatterplot/img/example-1.png "categoricalscatterplot")

### Example 2: Customizing CategoricalScatterplot plots
```
% A minimalistic categorical scatter plot
CategoricalScatterplot(MPG,cellstr(Origin),'WhiskerLine',false,'BoxColor',[0.8471 0.8627 0.8392],'WhiskerColor',[0.8235 0.7412 0.0392])
```
Categorical Scatter plots can easily be customized using the input arguments to the function. Here is a minimalistic version of the above plot:
![](https://github.com/AbstractGeek/CategoricalScatterplot/img/example-2.png "categoricalscatterplot")


### Syntax
#### Required Inputs:
X - Input Data - vector or a matrix (Group required if X is a vector)

#### Optional Inputs:
Group - Grouping variables - vector
'Color' - nx3 matrix for n groups or 1x3 vector for all groups or a
character color ('k', 'r', etc) for all groups
'Labels' - A cell string containing labels for all the groups

####Plot parameters (to tweak style):
*Scatter parameters:*
'binWidth': Sets the bin width that is used to stagger points along the
 xaxis in the scatter plot. Smaller binwidths imply that the points will
 be close to the central line. Larger values will make the points
 distributed all over the box.
'binWidthRatio': Easier way of setting binwidth. It calulates the
 binwidth automatically based on the value range of points (Y axis range).
'spreadWidth': Sets the extent of the point spread on the x axis.
'boxWidth': Sets the width of the boxes.
*Plotting styles:*
'Marker': Marker type - char
'MarkerSize': Marker size - num
'FillMarker': Logical (true or false / 0 or 1)
'BoxColor': char ('r') or rgb vector ([0 0 1])
'BoxEdgeColor': char ('r') or rgb vector ([0 0 1])
'MedianColor': char ('r') or rgb vector ([0 0 1])
'WhiskerColor': char ('r') or rgb vector ([0 0 1])
'BoxAlpha': Transparency of the box. num (0 to 1)
'BoxLineStyle': char ('-')
'MedianLineStyle': char ('-')
'WhiskerLineStyle': char ('-')
'BoxLineWidth': num (2.0)
'MedianLineWidth': num (2.0)
'WhiskerLineWidth': num (2.0)
