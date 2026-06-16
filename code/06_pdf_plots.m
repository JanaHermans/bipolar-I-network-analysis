% This script plots the probability density function (PDF) for each
% anxiety and mania subgroup

%% CAS clusters
load('CAS_RUN10_OCT_22_env.mat')

% get some dimensions for subplots based on groupsize
dim1 = length(unique(groups));
dim2 = ceil(sqrt(dim1));

figure; hold on;
for i = 1:dim1
un = unique(groups);
r = un(i,1);
x = squeeze(median(Adj(:,:,(groups == r)),3));
x = x(:);
x = x(~isnan(x));
x = x(x>0);
x = x(x<1);
[f,xi] = ksdensity(x);
if i == 1
x0 = xi;
else i == 2
x1 = xi;
end
area(xi,f,'FaceAlpha', 0.5);
title(strcat("Probability Density Function for Anxiety subgroups"), 'FontSize', 16);
end
hold off
legend(strcat("Cluster ",num2str(un(1,1))),strcat("Cluster ",num2str(un(2,1))))

%% YMRS clusters
load('YMRS_RUN10_OCT_22_env.mat')

% get some dimensions for subplots based on groupsize
dim1 = length(unique(groups));
dim2 = ceil(sqrt(dim1));

figure; hold on;
for i = 1:dim1
    un = unique(groups);
    r = un(i,1);
    x = squeeze(median(Adj(:,:,(groups == r)),3));
    x = x(:);
    x = x(~isnan(x));
    x = x(x>0);
    x = x(x<1);
    [f,xi] = ksdensity(x);
    if i == 1
        x0 = xi;
    elseif i == 2
        x1 = xi;
    elseif i == 3
        x2 = xi;
    else
        x3 = xi;
    end
    area(xi,f,'FaceAlpha', 0.5);
    title(strcat("Probability Density Function for Mania subgroups"), 'FontSize', 16);
end
hold off
legend(strcat("Cluster ",num2str(un(1,1))),strcat("Cluster ",num2str(un(2,1))),strcat("Cluster ",num2str(un(3,1))),strcat("Cluster ",num2str(un(4,1))))


