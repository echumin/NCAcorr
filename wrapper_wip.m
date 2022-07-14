%% Network Contingency Analysis - Correlation

% INPUTS
%  ts      - 3D Time by Node by Subject matrix
%  net     - System affiliation vector (e.g. yeo7 labels)
%  bins    - 2D binsets by number of bins matrix of ones and zeros, where
%            each row is a set of bins: 1=included 0=excluded root sum squared
%            (RSS) ordered data bin.
%            [1 1 1 1] -> include all data (equivalent to pearson FC)
%            [1 0 0 0] -> include only time points from lowest 25% RSS 
%            (first of four bins). 
%  scores  - 2D number of subjects by number of scores matrix for
%            correlation with RSS component data.
%  edgethr - [0.01] Edge level initial significance threshold.
%  edgextnt- [] Extent threshold. Default=no threshold. 
%            For values 0 to <1 : minimum fraction of edges in block that
%            must be significant for the block to be considered.
%            For values >1 minimum number of edges. 
%  nl_mod  - [1] 1 - null model scrambles the FC-cogscore pairings by
%                         permuting the scores over the sample.
%                2 - null model scrambels the block structure by
%                         permuting the block index. This is done post null
%                         correlation, so all its doing is randomly
%                         redistributing the significance.
%                3 - [overkill?] null model scrambles the unique
%                         EDGES independently, destroying the block
%                         structure. Additionally, permutation is done
%                         before computing correlation, independently for
%                         each subject, so any intersubject relationship is
%                         destroyed.
%  nPerm   - [5000] Number of score permutations for block level
%            significance.
%  multcm  - ['bonferroni' or 'fdr'] multiple comparison correction
%            strategy for block level significance.
%  corrP   - [0.05] multiple comparisons adjusted threshold.
%  group   - provide a vector of group assignments. (give all ones for full
%            sample for now)

% get Functional Connectivity (FC) components for an input number of bins
% from RSS ordered edge time series. 
% Time series are ordered by the input network affilications

%% Step 0: Add the NCAcorr package to path
    % I havent moved the wrapper, so working directory and its subdirectories
    % are added to path.
addpath(genpath(pwd))

%% Step1: Generate RSS based FC components for input bin sets.
    % NOTE that equal size binning is enforced. This means the remainder
    % that can't be disributed into all bins is divided by 2 and cropped
    % from the time series start and end.
[FCcomponents, block_idx] = fcn_binRSS(ts,net,bins);

%% Step2: Run NCAcorr
    % NOTE: Spearman correlation is hard coded.
       
% Example: ~FC and Pentile (20%) bins. 
bins = [1,1,1,1,1;
        1,0,0,0,0;
        0,1,0,0,0;
        0,0,1,0,0;
        0,0,0,1,0;
        0,0,0,0,1];
% ncc_out = fcn_netcontcorr(FCcomponents,bins,block_idx,scores,edgethr,edgextnt,nl_mod,nPerm,   multcm   ,corrP)    
  ncc_out = fcn_netcontcorr(FCcomponents,bins,block_idx,scores, .01,       [],     1,   5000,'bonferroni',.05)

% Network and score label examples:
netLabels = {'VIS','SOM','DAN','VAN','LIM','FRP','DMN'};
scoreLabels = {'AttnProcSpeed','CogFunc','Memory'};
% group -> subject group affilication vector (if applicable)

%% Step3: Create summary figures.
ncc_summary(ncc_out,scoreLabels,netLabels,group)
% whole sample variant:
ncc_summary(ncc_out,scoreLabels,netLabels,ones(size(ts,3)))




