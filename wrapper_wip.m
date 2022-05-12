%% Network Contingency Analysis - Correlation

% INPUTS
%   ts     - 3D Time by Node by Subject matrix
%   net    - System affiliation vector (e.g. yeo7 labels)
%   bins   - 2D binsets by number of bins matrix of ones and zeros, where
%            each row is a set of bins: 1=included 0=excluded root sum squared
%            (RSS) ordered data bin.
%            [1 1 1 1] -> include all data (equivalent to pearson FC)
%            [1 0 0 0] -> include only time points from lowest 25% RSS 
%            (first of four bins). 
%   scores - 2D number of subjects by number of scores matrix for
%            correlation with RSS component data.
%   thr    - [0.01] Edge level initial significance threshold.
%   nPerm  - [10000] Number of score permutations for block level
%            significance.
%   multcm - ['bonferroni' or 'fdr'] multiple comparison correction
%            strategy for block level significance.
%   corrP  - [0.05] multiple comparisons adjusted threshold.

% get Functional Connectivity (FC) components for an input number of bins
% from RSS ordered edge time series. 
% Time series are ordered by the input network affilications

addpath(genpath(pwd))

%% Step1: Generate RSS based FC components for input bin sets.
[FCcomponents, block_idx] = fcn_binRSS(ts,net,bins);

%% Step2: Run NCAcorr
ncc_out = fcn_netcontcorr(FCcomponents,bins,block_idx,scores,.01,25000,'bonferroni',.05)

% Network and score labels
netLabels = {'VIS','SOM','DAN','VAN','LIM','FRP','DMN'};
scoreLabels = {'AttnProcSpeed','CogFunc','Memory'};
% group -> subject group affilications (if applicable)

%% Step3: Create summary figures.
ncc_summary(ncc_out,scoreLabels,netLabels,group)



