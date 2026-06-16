% This script loads the connectivity matrices, applies a threshold to retain the
% top 10 percent of connections, and calculates global and regional graph
% measures for each subject.

function [] = wrapper_YMRS_function(name)
% Load group labels and set variables
cd '/home/'
cluster_list = readtable('newcom_label_ymrs_cas_mat_FINAL.csv');
subjects = cluster_list.('src_subject_id');
cluster_list = rmmissing(cluster_list);
% Load group labels 
ymrs_c = cluster_list.('YMRS_cluster');
groups = ymrs_c(:,:);
% Import matrices into 3D structure
cd '/home/'
for i = 1:length(groups)
    subID = string(subjects(i,1));
    cd(strcat('/home/sub-',subID))
    file = strcat('Deconfounded_cor_mat_sub-', subID, '.txt');
    matrix = readtable(file,'ReadRowNames',true, 'HeaderLines', 1);
    matrix = table2array(matrix);
    Adj(:,:,i) = matrix;
    cd '/home/'
end
% Set rest of variables
nRand = 10;
percentage = 0.1;
wResult.nRand = nRand;
regionLabels = num2cell([1:size(Adj,1)]);
wResult.group = groups;
wResult.regionLabels = regionLabels;
PlotGlobal = 1;
PlotLocal = 0;
PlotMatrices = 1;
for i = 1:size(Adj,3); filelist{i} = strcat('file_',num2str(i));end;
cd '/home/'
% Weighted networks
nSubjects = size(Adj,3);
n=size(Adj,1);
Adj(find(Adj<0))=0;
Adj(1:n+1:n*n)=1;
Adj_2 = Adj;
h = waitbar(0,'Running analysis on weighted networks');
for i = 1:nSubjects
    waitbar(i/nSubjects);
    A = squeeze(Adj(:,:,i));
    A = (A+A')./2.;
    G = graph(A);
    T = minspantree(G,'Method','sparse');
    MST = adjacency(T);
    MST = full(MST);
    TOP = threshold_proportional(A,percentage);
    TOP = double(TOP>0);
    COM = MST + TOP;
    COM = double(COM>0);
    A = times(COM,A);
    Adj_2(:,:,i) = A;    
    % create non-normalized output
    wResult.deg(i,:) = degrees_und(A); %degree
    wResult.dens(i,:) = density_und(A); %density
    wResult.cpl(i,:) = charpath(A); %characteristic path length
    wResult.trans(i,:) = transitivity_wu(A); %transitivity
    wResult.strength(i,:) = strengths_und(A); %strength
    wResult.assor(i,:) = assortativity_wei(A,0); %assortativity
    [M, Q] = modularityConsensusFun(A,1,nRand); %optimized clustering coefficient
    wResult.M(i,:) = M;
    wResult.clustcoeff(i,:) = mean(clustering_coef_bu(A));
    wResult.part(i,:) = participation_coef(A,M);    
    for iR = 1:nRand
        R = randmio_und_connected(A, 10);% create a random matrix from original
        cpl_random(iR,:) = charpath(distance_wei(R));
        clust_random(iR,:) = mean(clustering_coef_wu(R));
        trans_random(iR,:) = transitivity_bu(R);
    end    
    CplRand = mean(cpl_random(find(not(isinf(cpl_random)))));
    ClustRand = mean(clust_random,1);
    TransRand = mean(trans_random,1);    
    wResult.Norm.cpl(i,:) = wResult.cpl(i,:)/CplRand;
    wResult.Norm.clustcoeff(i,:) = wResult.clustcoeff(i,:)/ClustRand;
    wResult.Norm.trans(i,:) = wResult.trans(i,:)/TransRand;    
    [wResult.mask(i,:), wResult.net(:,:,i)] = getTop(wResult.deg(i,:),A,percentage);
end
close(h)
%% Mean connectivity
for i = 1:nSubjects
    A = squeeze(Adj(:,:,i));
    A = A(A>0);
    wResult.meancon(i,:) = mean(A(:)); %calculate the mean connectivity
end
wResult.subjects = subjects;
wResult.clusterlabels = cluster_list.('YMRS_cluster');
cd '/home/'
save([name '_env'])
save([name '_structure.mat'], '-struct', 'wResult')
end
