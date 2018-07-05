function options = lsmtb_options
% retrieve lesion symptom mapping toolbox options. this code should be
% refactored and this function made obsolete. One day soon...
options.lesionvolumecorrection = {'Regress on Behavior','Regress on Lesion','Regress on Both','DTLVC','None'};
%options.old_hypodirection = {'One-tailed (positive)','One-tailed (negative)','Two-tailed'}; % terms used before 2/7/18
%options.hypodirection = {'High scores are bad','High scores are good','Two-tailed'};
options.hypodirection = {'High scores are bad','High scores are good'}; % Two-tailed removed for first release since it wasn't implemented
