function [] = categoricalscatterplot(X, Group, varargin)
% function [] = categoricalscatterplot()
%
%
%
% Dinesh Natesan
% Last Modified: 16th Aug 2016

%% Parse Inputs
p = inputParser;
% Main arguments
addRequired(p, 'X', @(X) ismatrix(X));
addRequired(p, 'Group', @(X) ismatrix(X) || iscellstr(X));
addOptional(p, 'Color', 'k', @(X) all(size(X) == [length(unique(Group)), 3]) ||...
    all(size(X) == [1, 3]) || (ischar(X) && length(X)==1));
addOptional(p, 'Labels', false, @(X) iscellstr(X));

% Scatter parameters
addParameter(p, 'binWidth', false, @(x) isnumeric(x));
addParameter(p, 'binWidthRatio', 0.05, @(x) isnumeric(x));
addParameter(p, 'spreadWidth', 0.6, @(x) isnumeric(x));
addParameter(p, 'boxWidth', 0.6, @(x) isnumeric(x));
% Plotting styles
addParameter(p, 'Marker', 'o', @(X) (ischar(X) && length(X)==1));
addParameter(p, 'MarkerSize', 30, @(x) isnumeric(x));
addParameter(p, 'FillMarker', true, @(x) islogical(x));
addParameter(p, 'BoxColor', [0.8471    0.8627    0.8392], @(X) ...
    all(size(X) == [1, 3]) || (ischar(X) && length(X)==1));
addParameter(p, 'BoxEdgeColor', 'none', @(X) ...
    all(size(X) == [1, 3]) || (ischar(X) && (length(X)==1 || length(X)==4)));
addParameter(p, 'MedianColor', 'r', @(X) ...
    all(size(X) == [1, 3]) || (ischar(X) && length(X)==1));
addParameter(p, 'WhiskerColor', [0.8235    0.7412    0.0392], @(X) ...
    all(size(X) == [1, 3]) || (ischar(X) && length(X)==1));
addParameter(p, 'BoxAlpha', 0.50, @(x) isnumeric(x));
addParameter(p, 'BoxLineStyle', '-', ...
    @(X) ischar(X) && (length(X)==1 || length(X)==2));
addParameter(p, 'MedianLineStyle', '-', ...
    @(X) ischar(X) && (length(X)==1 || length(X)==2));
addParameter(p, 'WhiskerLineStyle', '-', ...
    @(X) ischar(X) && (length(X)==1 || length(X)==2));
addParameter(p, 'BoxLineWidth', 1.0, @(x) isnumeric(x));
addParameter(p, 'MedianLineWidth', 1.0, @(x) isnumeric(x));
addParameter(p, 'WhiskerLineWidth', 1.5, @(x) isnumeric(x));
% Parse inputs and unpack structure
parse(p, X, Group, varargin{:});
parsed = p.Results;

%% Convert the groups into a cell array
groups = unique(Group);
data = cell(length(groups), 2);
new_data = cell(length(groups), 2);
Xlen = 0;
Xmax = -999;
Xmin = 999;

for i = 1:length(groups)
    data{i,1} = X(Group == groups(i));
    data{i,2} = Group(Group == groups(i));
    if (Xlen < length(data{i,1}))
        Xlen = length(data{i,1});
    end
    
    if (Xmin > floor(min(data{i,1})))
        Xmin = floor(min(data{i,1}));
    end
    
    if (Xmax < ceil(max(data{i,1})))
        Xmax = ceil(max(data{i,1}));
    end
    
end

% Get binWidth ratio is only ratio is provided
if (islogical(parsed.binWidth) && ~parsed.binWidth)
    binWidth = parsed.binWidthRatio * round(Xmax - Xmin);
else
    binWidth = parsed.binWidth;
end

%% Discretize points in a group
for i = 1:length(groups)
    Xtemp = data{i,1};
    Ytemp = data{i,2};    
    [counts,~,bins] = histcounts(Xtemp, 'BinWidth', binWidth);
    inds = find(counts~=0);
    counts = counts(inds);
    
    for j=1:length(inds)
        width = parsed.spreadWidth * (1-exp(-0.1 * (counts(j)-1)));
        xpoints = linspace(-width/2, width/2, counts(j)) + i;
        Ytemp(bins==inds(j)) = xpoints;
    end
    
    new_data{i,1} = Xtemp;
    new_data{i,2} = Ytemp;
    
end

%% Plot the data beautifully
boxWidth = parsed.boxWidth;
hold on;

for i = 1:length(groups)    
    imp_quantiles = quantile(new_data{i,1}, [0.25, 0.5, 0.75]);
    IQR = imp_quantiles(3) - imp_quantiles(1);
    whisker = [imp_quantiles(1) - 1.5 * IQR, ...
        imp_quantiles(3) + 1.5 * IQR];
    
    % Draw box
    patch([i-boxWidth/2, i-boxWidth/2, i+boxWidth/2, i+boxWidth/2]',...
        [imp_quantiles(3), imp_quantiles(1), imp_quantiles(1), imp_quantiles(3)]',...
        parsed.BoxColor, 'FaceAlpha', parsed.BoxAlpha,...
        'EdgeColor', parsed.BoxEdgeColor, ...
        'LineStyle', parsed.BoxLineStyle, 'LineWidth', parsed.BoxLineWidth);
    
    % Draw points
    if parsed.FillMarker
        if (ischar(parsed.Color)||size(parsed.Color,1))
            scatter(new_data{i,2}, new_data{i,1}, parsed.MarkerSize, 'filled',...
                parsed.Color);
        else
            scatter(new_data{i,2}, new_data{i,1}, parsed.MarkerSize, 'filled',...
                parsed.Color(i,:));
        end
    else
        if (ischar(parsed.Color)||size(parsed.Color,1))
            scatter(new_data{i,2}, new_data{i,1}, parsed.MarkerSize, 'filled',...
                parsed.Color);
        else
            scatter(new_data{i,2}, new_data{i,1}, parsed.MarkerSize, 'filled',...
                parsed.Color(i,:));
        end
    end
    
    % Draw median
    plot([i-boxWidth/2, i+boxWidth/2], ...
        [imp_quantiles(2), imp_quantiles(2)], ...
        'LineStyle', parsed.MedianLineStyle, ...
        'Color', parsed.MedianColor,...
        'LineWidth', parsed.MedianLineWidth);
    
    % Draw Q + 1.5 IQR
    temp = sortrows([whisker(2) - data{i,1}, (1:length(data{i,1}))'], 1);
    closest_point = temp(find(temp(:,1) >=0, 1, 'first'),2);
    plot([i-boxWidth/5, i+boxWidth/5], ...
        [new_data{i,1}(closest_point), new_data{i,1}(closest_point)], ...
        'LineStyle', parsed.WhiskerLineStyle, ...
        'Color', parsed.WhiskerColor,...
        'LineWidth', parsed.WhiskerLineWidth);
    
    % Draw Q - 1.5 IQR
    temp = sortrows([data{i,1} - whisker(1), (1:length(data{i,1}))'], 1);
    closest_point = temp(find(temp(:,1) >=0, 1, 'first'),2);
    plot([i-boxWidth/5, i+boxWidth/5], ...
        [new_data{i,1}(closest_point), new_data{i,1}(closest_point)], ...
        'LineStyle', parsed.WhiskerLineStyle, ...
        'Color', parsed.WhiskerColor,...
        'LineWidth', parsed.WhiskerLineWidth);
    
end

ax = gca;
ax.XTick = 1:length(groups);
if (islogical(parsed.Labels) && ~parsed.Labels)
    ax.XTickLabel = groups;
else
    ax.XTickLabel = parsed.Labels;
end

end