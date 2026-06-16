% This script plots the mean connectivity matrices (thresholded and unthresholded)
% for each anxiety and mania subgroup.

%% CAS matrix plots

load('CAS_RUN10_OCT_22_env.mat')

% get some dimensions for subplots based on groupsize
dim1 = 2;
dim2 = 1;

if PlotMatrices == 1
    figure; hold on;
    for i = 1:dim1
        un = unique(groups);
        p = un(i,1);
        mat = squeeze(median(Adj(:,:,(groups == p)),3));
        
        subplot(dim1/dim2,ceil(dim1/dim2),i);
        imagesc(mat);            % Create a colored plot of the matrix values
        title(strcat("Unthresholded connectivity matrix for anxiety subgroup: ",num2str(p)), 'FontSize', 10);
        colormap(parula);  % Change the colormap to gray (so higher values are
        %#   black and lower values are white)
        c = colorbar; ylabel(c,'median connectivity score ')
        caxis([0.0 1.0])
        
        %         textStrings = num2str(mat(:),'%0.2f');  %# Create strings from the matrix values
        %         textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
        %         [x,y] = meshgrid(1:size(mat,1));   %# Create x and y coordinates for the strings
        %         hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
        %             'HorizontalAlignment','center');
        %         midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
        %         textColors = repmat(mat(:) > midValue,1,3);  %# Choose white or black for the
        %         %#   text color of the strings so
        %         %#   they can be easily seen over
        %         %#   the background color
        %         set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
        %
        %         set(gca,'XTick',1:size(mat,1),...                         %# Change the axes tick marks
        %             'XTickLabel',regionLabels,...  %#   and tick labels
        %             'YTick',1:size(mat,1),...
        %             'YTickLabel',regionLabels,...
        %             'TickLength',[0 0]);
        
        
    end
    hold off
    
end

if PlotMatrices == 1
    figure; hold on;
    for i = 1:dim1
        un = unique(groups);
        p = un(i,1);
        mat = squeeze(median(Adj_2(:,:,(groups == p)),3));
        
        subplot(dim1/dim2,ceil(dim1/dim2),i);
        imagesc(mat);            % Create a colored plot of the matrix values
        title(strcat("Thresholded connectivity matrix for anxiety subgroup: ",num2str(p)), 'FontSize', 10);
        colormap(parula);  % Change the colormap to gray (so higher values are
        %#   black and lower values are white)
        c = colorbar; ylabel(c,'median connectivity score ')
        caxis([0.0 1.0])
        
        %         textStrings = num2str(mat(:),'%0.2f');  %# Create strings from the matrix values
        %         textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
        %         [x,y] = meshgrid(1:size(mat,1));   %# Create x and y coordinates for the strings
        %         hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
        %             'HorizontalAlignment','center');
        %         midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
        %         textColors = repmat(mat(:) > midValue,1,3);  %# Choose white or black for the
        %         %#   text color of the strings so
        %         %#   they can be easily seen over
        %         %#   the background color
        %         set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
        %
        %         set(gca,'XTick',1:size(mat,1),...                         %# Change the axes tick marks
        %             'XTickLabel',regionLabels,...  %#   and tick labels
        %             'YTick',1:size(mat,1),...
        %             'YTickLabel',regionLabels,...
        %             'TickLength',[0 0]);
        
        
    end
    hold off
    
end
%% YMRS matrix plots

load('YMRS_RUN10_OCT_22_env.mat')

% get some dimensions for subplots based on groupsize
dim1 = length(unique(groups));
dim2 = ceil(sqrt(dim1));

if PlotMatrices == 1
    figure; hold on;
    for i = 1:dim1
        un = unique(groups);
        p = un(i,1);
        mat = squeeze(median(Adj(:,:,(groups == p)),3));
        
        subplot(dim1/dim2,ceil(dim1/dim2),i);
        imagesc(mat);            % Create a colored plot of the matrix values
        title(strcat("Unthresholded connectivity matrix for mania subgroup: ",num2str(p)), 'FontSize', 10);
        colormap(parula);  % Change the colormap to gray (so higher values are
        %#   black and lower values are white)
        c = colorbar; ylabel(c,'median connectivity score ')
        caxis([0.0 1.0])
        
        %         textStrings = num2str(mat(:),'%0.2f');  %# Create strings from the matrix values
        %         textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
        %         [x,y] = meshgrid(1:size(mat,1));   %# Create x and y coordinates for the strings
        %         hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
        %             'HorizontalAlignment','center');
        %         midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
        %         textColors = repmat(mat(:) > midValue,1,3);  %# Choose white or black for the
        %         %#   text color of the strings so
        %         %#   they can be easily seen over
        %         %#   the background color
        %         set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
        %
        %         set(gca,'XTick',1:size(mat,1),...                         %# Change the axes tick marks
        %             'XTickLabel',regionLabels,...  %#   and tick labels
        %             'YTick',1:size(mat,1),...
        %             'YTickLabel',regionLabels,...
        %             'TickLength',[0 0]);
        
        
    end
    hold off
    
end

if PlotMatrices == 1
    figure; hold on;
    for i = 1:dim1
        un = unique(groups);
        p = un(i,1);
        mat = squeeze(median(Adj_2(:,:,(groups == p)),3));
        
        subplot(dim1/dim2,ceil(dim1/dim2),i);
        imagesc(mat);            % Create a colored plot of the matrix values
        title(strcat("Thresholded connectivity matrix for mania subgroup: ",num2str(p)), 'FontSize', 10);
        colormap(parula);  % Change the colormap to gray (so higher values are
        %#   black and lower values are white)
        c = colorbar; ylabel(c,'median connectivity score ')
        caxis([0.0 1.0])
        
        %         textStrings = num2str(mat(:),'%0.2f');  %# Create strings from the matrix values
        %         textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
        %         [x,y] = meshgrid(1:size(mat,1));   %# Create x and y coordinates for the strings
        %         hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
        %             'HorizontalAlignment','center');
        %         midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
        %         textColors = repmat(mat(:) > midValue,1,3);  %# Choose white or black for the
        %         %#   text color of the strings so
        %         %#   they can be easily seen over
        %         %#   the background color
        %         set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
        %
        %         set(gca,'XTick',1:size(mat,1),...                         %# Change the axes tick marks
        %             'XTickLabel',regionLabels,...  %#   and tick labels
        %             'YTick',1:size(mat,1),...
        %             'YTickLabel',regionLabels,...
        %             'TickLength',[0 0]);
        
        
    end
    hold off
    
end
