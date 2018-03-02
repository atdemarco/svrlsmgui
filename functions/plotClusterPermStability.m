function out = plotClusterPermStability(varargin)
% This function will plot the stability of cluster size for an svrlsmgui
% run given an analysis clusterwisedir at each multiple of interval assess_interval
% and plot the top ntopclusters clusters in the actual observed data
% ... input parameters -- 1 analysis rootdir, 2 assess interval
% AD 8/4/17
% added getframe on 9/1/17
% added references to specific figures (invis) for plot commands 1/20/18
% as of 2/14/18, the first input is the actual cluster file.
if nargin > 0
    %clusterwisedir = varargin{1};
    parms = varargin{1};
else % default for testing...
    %error('One input is required pointing toward the cluster dir');
    error('First input should be parms struct for analysis to assess.');
end
if nargin > 1
    assess_interval = varargin{2};
else
    assess_interval = 100;
end
if nargin > 2
    ntopclusters = varargin{3};
else
    ntopclusters = 5; % plot top 5 clusters...
end

% get requested critical value from directory name 
% [~,thisdir]=fileparts(clusterwisedir); 
% words = strsplit(thisdir);
% pvalstr = ['.' words{2}(2:end)];
% pval = str2num(pvalstr);
pval = parms.clusterwise_p;
pvalstr = num2str(pval);
critical_percentile = round(100*(1-pval));

% read in all max cluster size list
%clustdata = load(fullfile(clusterwisedir,'Largest clusters.mat'));
clustdata = load(parms.files_created.largest_clusters);
clustdata = clustdata.all_max_cluster_sizes;

nclusterperms = numel(clustdata);

vals = []; 
assess_indices = [(assess_interval : assess_interval : nclusterperms)] ;
for a = assess_indices
    relevant_vals = clustdata(1:a); 
    vals(end+1) = prctile(relevant_vals,critical_percentile);
end

clusterfig = figure('visible','off');
a = axes(clusterfig);
colororder = repmat('rcmgby',1,5);
plot(assess_indices,vals,'k.-','parent',a)
xlims = get(a,'xlim');

% Now show observed clusters..
%empirically_observed_clusters = readtable(fullfile(clusterwisedir,'Table of clusters.txt'));
empirically_observed_clusters = readtable(parms.files_created.clustertable);
legendvals = {['Critical value ' pvalstr]}; % we'll add onto this with cluster identities.
for c = 1 : ntopclusters
    if size(empirically_observed_clusters,1) >= c % make sure it exists.
        curclustersize = empirically_observed_clusters.nvox(c);
        line(xlims,[curclustersize curclustersize],'Color',colororder(c),'LineStyle','--','parent',a);
        legendvals{end+1} = ['Cluster ' num2str(c) ' (' num2str(curclustersize) ')'];
    end
    
end
legend(legendvals)
title(['Stability of crit clust vol (P = ' pvalstr ') by perm, N = ' num2str(nclusterperms)],'parent',a)
ylabel('Cluster volume','parent',a)
xlabel('Permutation number','parent',a)

F = getframe(clusterfig); % capture whole figure.
close(clusterfig);
out = F.cdata;

