function [ncc_out] = fcn_netcontcorr(FCcomponents,binset,block_idx,cogscores,edge_thr,edge_xtnt,nl_mod,nPerm,multccorr,mc_thr)
%%                   Network Contingency Correlation                     %%
% INPUTS
%   FCcomponents -> Data subset components based on binSets 
%   binSet       -> 2D binsets by number of bins matrix of ones and zeros, where
%                   each row is a set of bins: 1=included 0=excluded root sum squared
%                   (RSS) ordered data bin.
%                   [1 1 1 1] -> include all data (equivalent to pearson FC)
%                   [1 0 0 0] -> include only time points from lowest 25% RSS 
%                   (first of four bins). 
%   block_idx    -> Node system affiliations that are consistent with
%                   reordered time series.
%   cogscores    -> 2D number of subjects by number of scores matrix for
%                   correlation with RSS component data.
%   edge_thr     -> [0.01] Edge level initial significance threshold.
%   edge_xtnt    -> [] Extent threshold. Default=no threshold. 
%                   For values 0 to <1 : minimum fraction of edges in block that
%                   must be significant for the block to be considered.
%                   For values >1 minimum number of edges. 
%   nl_mod       -> [1] 1 - null model scrambles the FC-cogscore pairings by
%                         permuting the scores over the sample.
%                       2 - null model scrambels the block structure by
%                         permuting the block index. This is done post null
%                         correlation, so all its doing is randomly
%                         redistributing the significance.
%                       3 - [overkill?] null model scrambles the unique
%                         EDGES independently, destroying the block
%                         structure. Additionally, permutation is done
%                         before computing correlation, independently for
%                         each subject, so any intersubject relationship is
%                         destroyed.
%  nPerm         -> [10000] Number of score permutations for block level
%                   significance.
%  multccorr     -> ['bonferroni' or 'fdr'] multiple comparison correction
%                   strategy for block level significance.
%  mc_thr        -> [0.05] multiple comparisons adjusted threshold.

% OUTPUTS
%   ncc_out -> a structure of input/output variables needed for the summary
%   code.
%
% Evgeny Jenya Chumin, 2023, Indiana University
%%
narginchk(4,10)
% edges, subjects, bins
[K, S, B] = size(FCcomponents);

if B ~= size(binset,1)
    error('Number of FC component bins does not match nunber of binsets');
end
ncc_out.nBins = size(binset,2);
ncc_out.binset = binset;
N = roots([-.5,.5,K]);
N = int16(N(N>0));
if N~=length(block_idx)
    error('Number of nodes in block_idx is different than in FCcomponents');
end
inds = find(triu(ones(N),1));
if ~exist('edge_thr','var')
    edge_thr = 0.01;
end
if ~exist('edge_xtnt','var')
    edge_xtnt = [];
    disp('No edge extent threshold enforced.')
elseif edge_xtnt <= 1
    disp('Edge extent threshold set as franction of edges per block')
elseif edge_xtnt > 1
    disp('Edge extent threshold set as absolute number of edges.')
elseif isempty(edge_xtnt)
    disp('No edge extent threshold enforced.')
end
if ~exist('nl_mod','var')
    nl_mod = 1;
elseif nl_mod == 1
    disp('Score permuting null model will be used.')
elseif nl_mod == 2
    disp('Block permuting null model will be used.')
elseif nl_mod == 3
    disp('Edge permuting null model will be used')
else
    error('Unrecognized null model argument')
end 
if ~exist('nPerm','var')
    nPerm = 5000;
end
if ~exist('mc_thr','var')
    mc_thr = 0.05;
end
nSys = max(block_idx);
if size(cogscores,1) ~= S
    error('Number of subjects between FCcomponents and cogscores inconsistent');
end
ncc_out.block_idx = block_idx;
ncc_out.cogscores = cogscores;
ncc_out.FCcomponents = FCcomponents;

for bs=1:B % for every bin
    nFCc_vec_all = FCcomponents(:,:,bs)';
   
    for sc=1:size(cogscores,2) % for every cog score
        bscore = cogscores(:,sc);
        % correlate each edge with score (empirical)
        [BBcc,BBp] = corr(nFCc_vec_all,bscore,'rows','complete','type','s');
        
        temp = zeros(N); temp(inds) = squeeze(BBp); BBmat = temp+temp';
        % appply edge significance threshold (empirical)
        mask = (BBmat<edge_thr).*~eye(N);
        
        [ncc_out.component_idx{bs,sc},cpsz_emp] = get_components(mask);
        ncc_out.component_sizes{bs,sc} = cpsz_emp;
        
        % apply blockwise edge extent threshold
        mdx = module_density(mask,block_idx,1,'sum');
        if edge_xtnt<=1
            mdxfull = module_density(ones(N),block_idx,1,'sum');
            mdxfull = round(mdxfull.*edge_xtnt);
            mdx(mdx<mdxfull)=0;
        elseif edge_xtnt > 1
            mdx(mdx<=edge_xtnt)=0;
        end

        mask_nnz = nnz(mask);
        
        % random permutations
        mdsum = zeros(nSys); cpszr_max = zeros(nPerm,1);
        if nl_mod == 3
            bmat=zeros(N);
            count=0;
            for i=1:7
                for j=1:7
                    count=count+1;
                    bmat(block_idx==i,block_idx==j)=count;
                end
            end
            bmat_vec=bmat(inds);
        end
        
        switch nl_mod
            case 1
                parfor r=1:nPerm
                    [~,BBpr] = corr(nFCc_vec_all,bscore(randperm(S)),'rows','complete','type','s');
                    temp = zeros(N); temp(inds) = BBpr; BBmat = temp+temp';
                    maskr = (BBmat<edge_thr).*~eye(N);
                    [~,cpsz] = get_components(maskr);
                    cpszr_max(r,1) = max(cpsz);  
                    mdr = module_density(maskr,block_idx,1,'sum');
                    mdsum = mdsum+(mdx>mdr);
                    maskr_nnz(r,1) = nnz(maskr);
                end
            case 2
                parfor r=1:nPerm
                    % This gives the same outcome as null strategy 1 NBS,
                    % so its commented out here.
%                     [~,BBpr] = corr(nFCc_vec_all,bscore,'rows','complete','type','s');
%                     temp = zeros(N); temp(inds) = BBpr; BBmat = temp+temp';
%                     maskr = (BBmat<edge_thr).*~eye(N);
%                     [~,cpsz] = get_components(maskr);
                    mdr = module_density(mask,block_idx(randperm(N)),1,'sum');
                    mdsum = mdsum+(mdx>mdr);
                end
                cpszr_max = 'null';  
                maskr_nnz = 'null';
            case 3
                for r=1:nPerm
                    parfor s=1:S
                        [~,perm_idx] = sort(bmat_vec(randperm(K)));
                        perm_nFCc_vec_all(s,:)=nFCc_vec_all(s,perm_idx);
                    end   
                    [~,BBpr] = corr(perm_nFCc_vec_all,bscore,'rows','complete','type','s');
                    temp = zeros(N); temp(inds) = BBpr; BBmat = temp+temp';
                    maskr = (BBmat<edge_thr).*~eye(N);
                    [~,cpsz] = get_components(maskr);
                    cpszr_max(r,1) = max(cpsz);  
                    mdr = module_density(maskr,block_idx,1,'sum');
                    mdsum = mdsum+(mdx>mdr);
                    maskr_nnz(r,1) = nnz(maskr);
                end
        end
        
        % p-vals
        mds_pval = 1-(mdsum./nPerm);
        if nl_mod ~=2
            cpsz_pval = 1-(sum((cpsz_emp>cpszr_max),1)/nPerm);
        else
            cpsz_pval = 'null';
        end
        
        % package
        
        ncc_out.modules_permP(:,:,bs,sc) = mds_pval;
        
        ncc_out.component_permP{bs,sc} = cpsz_pval;
        
        ncc_out.null_largest_components(:,bs,sc) = cpszr_max;
        
        ncc_out.edge_uncorrPthr_mask(:,:,bs,sc) = mask;
        
        ncc_out.num_uncorrPthr_edges(bs,sc) = mask_nnz;
        
        ncc_out.null_num_uncorrPthr_edges(:,bs,sc) = maskr_nnz;
        
        temp = zeros(N); temp(inds) = squeeze(BBp); BBpmat = temp+temp';
        ncc_out.BB_uncorrP_mat(:,:,bs,sc) = BBpmat;
        
        temp = zeros(N); temp(inds) = squeeze(BBcc); BBccmat = temp+temp';
        ncc_out.BB_corr_mat(:,:,bs,sc) = BBccmat;
        
        switch multccorr
            case 'bonferroni'
                ncc_out.corr_type = multccorr;
                zerosys = sum(mdx(logical(triu(ones(nSys))))==0);
                % only systems with nonzero uncorrected counts are tested
                ncc_out.corr_thr = mc_thr/(((nSys^2-nSys)/2+nSys)-zerosys);
                ncc_out.plot(bs,sc) = nnz(triu(mds_pval<=ncc_out.corr_thr));
            case 'fdr'
                ncc_out.corr_type = multccorr;
                mds_pval(mdx==0)=nan; % FDR will ignore nans in the computation
                ncc_out.corr_thr(bs,sc) = FDR(mds_pval(logical(triu(ones(nSys)))),mc_thr);
                ncc_out.plot(bs,sc) = nnz(triu(mds_pval<=ncc_out.corr_thr(bs,sc)));
            otherwise
                fprintf(2,'multccorr does not match existing options. Defaulting to Bonferoni.\n')
                ncc_out.corr_type = 'bonferroni';
                zerosys = sum(mdx(logical(triu(ones(nSys))))==0);
                % only systems with nonzero uncorrected counts are tested
                ncc_out.corr_thr = mc_thr/(((nSys^2-nSys)/2+nSys)-zerosys);
                ncc_out.plot(bs,sc) = nnz(triu(mds_pval<=ncc_out.corr_thr));
        end
        
        if bs==1 && sc==1
            disp(['binset ','cscore ','nsig'])
        end
        disp(['  ',num2str(bs),'      ',num2str(sc),'     ',num2str(ncc_out.plot(bs,sc))])
    end
end
