function [] = CategoricalScatterplot(X, varargin)
% function [] = categoricalscatterplot()
%
% A replacement for the traditional box and whisker plots provided in
% MATLAB (command boxplot). The categorical scatter plots additionally
% shows the data points, which is useful the visualize the underlying
% distribution (similar to violin plots).
%
% Required Inputs:
% X - Input Data - vector or a matrix (Group required if X is a vector)
%
% Optional Inputs:
% Group - Grouping variables - vector
% 'Color' - nx3 matrix for n groups or 1x3 vector for all groups
% 'Labels' - A cell string containing labels for all the groups
%
% Plot parameters (to tweak style):
% % Scatter parameters:
% 'binWidth': Sets the bin width that is used to stagger points along the
%  xaxis in the scatter plot. Smaller binwidths imply that the points will
%  be close to the central line. Larger values will make the points
%  distributed all over the box.
% 'binWidthRatio': Easier way of setting binwidth. It calulates the
%  binwidth automatically based on the value range of points (Y axis range).
% 'spreadWidth': Sets the extent of the point spread on the x axis.
% 'boxWidth': Sets the width of the boxes.
% % Plotting styles:
% 'Marker': Marker type - char
% 'MarkerSize': Marker size - num
% 'FillMarker': Logical (true or false / 0 or 1)
% 'BoxColor': char ('r') or rgb vector ([0 0 1])
% 'BoxEdgeColor': char ('r') or rgb vector ([0 0 1])
% 'MedianColor': char ('r') or rgb vector ([0 0 1])
% 'WhiskerColor': char ('r') or rgb vector ([0 0 1])
% 'BoxAlpha': Transparency of the box. num (0 to 1)
% 'BoxLineStyle': char ('-')
% 'MedianLineStyle': char ('-')
% 'WhiskerLineStyle': char ('-')
% 'BoxLineWidth': num (2.0)
% 'MedianLineWidth': num (2.0)
% 'WhiskerLineWidth': num (2.0)
%
% Dinesh Natesan (AbstractGeek)
% Last Modified: 11th Dec 2016

%% Parse Inputs
p = inputParser;
% Main arguments
addRequired(p, 'X', @(X) ismatrix(X));
addOptional(p, 'Group', [], @(X) isnumeric(X) || iscellstr(X));
addOptional(p, 'Color', [0 0 0], @(x) all(size(x) == size(X)) ||...
    all([size(X,1),3] == size(x)) || all(size(x) == [1, 3]));
addOptional(p, 'Labels', false, @(X) iscellstr(X));

% Scatter parameters
addParameter(p, 'binWidth', false, @(x) isnumeric(x));
addParameter(p, 'binWidthRatio', 0.05, @(x) isnumeric(x));
addParameter(p, 'spreadWidth', 0.6, @(x) isnumeric(x));
addParameter(p, 'boxWidth', 0.6, @(x) isnumeric(x));
% Plotting styles
addParameter(p, 'Marker', 'o', @(X) (ischar(X) && length(X)==1));
addParameter(p, 'MarkerSize', 25, @(x) isnumeric(x));
addParameter(p, 'FillMarker', true, @(x) islogical(x));
addParameter(p, 'WhiskerLine', true, @(x) islogical(x));
addParameter(p, 'BoxColor', [0.31, 0.31, 0.31], @(X) ...
    all(size(X) == [1, 3]) || (ischar(X) && length(X)==1));
addParameter(p, 'BoxEdgeColor', 'none', @(X) ...
    all(size(X) == [1, 3]) || (ischar(X) && (length(X)==1 || length(X)==4)));
addParameter(p, 'MedianColor', 'r', @(X) ...
    all(size(X) == [1, 3]) || (ischar(X) && length(X)==1));
addParameter(p, 'WhiskerColor', [0 0 0], @(X) ...
    all(size(X) == [1, 3]) || (ischar(X) && length(X)==1));
addParameter(p, 'BoxAlpha', 0.50, @(x) isnumeric(x));
addParameter(p, 'BoxLineStyle', '-', ...
    @(X) ischar(X) && (length(X)==1 || length(X)==2));
addParameter(p, 'MedianLineStyle', '-', ...
    @(X) ischar(X) && (length(X)==1 || length(X)==2));
addParameter(p, 'WhiskerLineStyle', '--', ...
    @(X) ischar(X) && (length(X)==1 || length(X)==2));
addParameter(p, 'BoxLineWidth', 1.0, @(x) isnumeric(x));
addParameter(p, 'MedianLineWidth', 1.0, @(x) isnumeric(x));
addParameter(p, 'WhiskerLineWidth', 1.0, @(x) isnumeric(x));
% Alternate Y axis
addParameter(p, 'YYaxis', false, @(X) ismatrix(X) || iscellstr(X) || ischar(X));
% Parse inputs and unpack structure
parse(p, X, varargin{:});
parsed = p.Results;
Group = parsed.Group;

% Handle color input
if all(size(parsed.Color)==[1,3])
    if isempty(Group)
        Color = repmat({parsed.Color}, size(X,1), size(X,2));
    else
        Color = repmat(parsed.Color, length(X), 1);
    end
else
    Color = parsed.Color;
end

% Handle fill input
if all(size(parsed.FillMarker) == [1,1])
    FillMarker = repmat(parsed.FillMarker,size(X,1),size(X,2));
else
    FillMarker = parsed.FillMarker;
end

if (size(X,2) == 1) || (size(X,1) == 1)
    %% Convert the groups into a cell array
    % Ensure group exists
    if isempty(Group)
        error('Group input required if X input is a vector');
    end
    
    if ischar(Group)
        Group = cellstr(Group);
    end
    
    % Make it into column vectors
    if (size(X,1) == 1)
        X = X';
    end
    
    % handle groups
    [group_names,~,group_ind] = unique(Group);
    groups = unique(group_ind);
    % Sort data into a cell array
    data = cell(length(groups), 2);
    new_data = cell(length(groups), 2);
    plot_vars = cell(length(groups), 2);
    Xmax = -999;
    Xmin = 999;
    
    for i = 1:length(groups)
        data{i,1} = X(group_ind == groups(i));
        data{i,2} = group_ind(group_ind == groups(i));
        plot_vars{i,1} = Color(group_ind == groups(i),:);
        plot_vars{i,2} = FillMarker(group_ind == groups(i),:);
        
        if (Xmin > floor(min(data{i,1})))
            Xmin = floor(min(data{i,1}));
        end
        
        if (Xmax < ceil(max(data{i,1})))
            Xmax = ceil(max(data{i,1}));
        end
        
    end
else
    %% Convert the matrix into cell array (after removing nans)
    group_names = Group;
    groups = 1:size(X,2);
    % Sort data into a cell array
    data = cell(length(groups), 2);
    new_data = cell(length(groups), 2);
    plot_vars = cell(length(groups), 2);
    Xmax = -999;
    Xmin = 999;
    
    for i=1:length(groups)
        data{i,1} = X(~isnan(X(:,i)),i);
        data{i,2} = i*ones(size(data{i,1}));
        plot_vars{i,1} = cell2mat(Color(~isnan(X(:,i)),i));
        plot_vars{i,2} = FillMarker(~isnan(X(:,i)),i);
        
        if (Xmin > floor(min(data{i,1})))
            Xmin = floor(min(data{i,1}));
        end
        
        if (Xmax < ceil(max(data{i,1})))
            Xmax = ceil(max(data{i,1}));
        end
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
    
    if parsed.YYaxis
        if (any(ismember(parsed.YYaxis,groups(i))))
            yyaxis right;
        else
            yyaxis left;
        end
    end
    
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
    scatter(new_data{i,2}(plot_vars{i,2}), new_data{i,1}(plot_vars{i,2}),...
        parsed.MarkerSize, plot_vars{i,1}(plot_vars{i,2},:), 'filled');
    scatter(new_data{i,2}(~plot_vars{i,2}), new_data{i,1}(~plot_vars{i,2}),...
        parsed.MarkerSize, plot_vars{i,1}(~plot_vars{i,2},:));
    
    % Draw median
    plot([i-boxWidth/2, i+boxWidth/2], ...
        [imp_quantiles(2), imp_quantiles(2)], ...
        'LineStyle', parsed.MedianLineStyle, ...
        'Color', parsed.MedianColor,...
        'LineWidth', parsed.MedianLineWidth,...
        'Marker','none');
    
    % Draw Q + 1.5 IQR
    if any(data{i,1}>whisker(2))
        % Outlier present
        plot([i-boxWidth/5, i+boxWidth/5], ...
            [whisker(2), whisker(2)], ...
            'LineStyle', '-', ...
            'Color', parsed.WhiskerColor,...
            'LineWidth', parsed.WhiskerLineWidth,...
            'Marker','none');
        
        if parsed.WhiskerLine
            % Draw top whiskers
            plot([i, i], [whisker(2), imp_quantiles(3)],...
                'LineStyle', parsed.WhiskerLineStyle, ...
                'Color', parsed.WhiskerColor,...
                'LineWidth', parsed.WhiskerLineWidth,...
                'Marker','none');
        end
        
    else
        % No outliers
        temp = sortrows([whisker(2) - data{i,1}, (1:length(data{i,1}))'], 1);
        closest_point = temp(find(temp(:,1) >=0, 1, 'first'),2);
        plot([i-boxWidth/5, i+boxWidth/5], ...
            [new_data{i,1}(closest_point), new_data{i,1}(closest_point)], ...
            'LineStyle', '-', ...
            'Color', parsed.WhiskerColor,...
            'LineWidth', parsed.WhiskerLineWidth,...
            'Marker','none');
        
        if parsed.WhiskerLine
            % Draw top whiskers
            plot([i, i], [new_data{i,1}(closest_point), imp_quantiles(3)],...
                'LineStyle', parsed.WhiskerLineStyle, ...
                'Color', parsed.WhiskerColor,...
                'LineWidth', parsed.WhiskerLineWidth,...
                'Marker','none');
        end
    end
    
    % Draw Q - 1.5 IQR
    if any(data{i,1}<whisker(1))
        % Outlier present
        plot([i-boxWidth/5, i+boxWidth/5], ...
            [whisker(1), whisker(1)], ...
            'LineStyle', '-', ...
            'Color', parsed.WhiskerColor,...
            'LineWidth', parsed.WhiskerLineWidth,...
            'Marker','none');
        
        if parsed.WhiskerLine
            % Draw top whiskers
            plot([i, i], [whisker(1), imp_quantiles(3)],...
                'LineStyle', parsed.WhiskerLineStyle, ...
                'Color', parsed.WhiskerColor,...
                'LineWidth', parsed.WhiskerLineWidth,...
                'Marker','none');
        end
        
    else
        % No outliers
        temp = sortrows([data{i,1} - whisker(1), (1:length(data{i,1}))'], 1);
        closest_point = temp(find(temp(:,1) >=0, 1, 'first'),2);
        plot([i-boxWidth/5, i+boxWidth/5], ...
            [new_data{i,1}(closest_point), new_data{i,1}(closest_point)], ...
            'LineStyle', '-', ...
            'Color', parsed.WhiskerColor,...
            'LineWidth', parsed.WhiskerLineWidth,...
            'Marker','none');
        
        if parsed.WhiskerLine
            % Draw bottom whiskers
            plot([i, i], [new_data{i,1}(closest_point), imp_quantiles(1)],...
                'LineStyle', parsed.WhiskerLineStyle, ...
                'Color', parsed.WhiskerColor,...
                'LineWidth', parsed.WhiskerLineWidth,...
                'Marker','none');
        end
    end
    
end

ax = gca;
ax.XTick = 1:length(groups);
if (islogical(parsed.Labels) && ~parsed.Labels)
    ax.XTickLabel = group_names;
else
    ax.XTickLabel = parsed.Labels;
end

end