function [ncc_out] = fcn_netcontcorr(FCcomponents,binset,block_idx,cogscores,edge_thr,nPerm,multccorr,mc_thr)
%%                   Network Contingency Correlation                     %%
% INPUTS
%   FCcomponents -> Data subset components based on binSets 
%   binSet - 2D binsets by number of bins matrix of ones and zeros, where
%            each row is a set of bins: 1=included 0=excluded root sum squared
%            (RSS) ordered data bin.
%            [1 1 1 1] -> include all data (equivalent to pearson FC)
%            [1 0 0 0] -> include only time points from lowest 25% RSS 
%            (first of four bins). 
%   block_idx    -> Node system affiliations that are consistent with
%                   reordered time series.


% OUTPUTS
%   ncc_out -> a structural input/output variables needed for the summary
%   code.


narginchk(4,8)
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
if ~exist('nPerm','var')
    nPerm = 10000;
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

for bs=1:B
    nFCc_vec_all = FCcomponents(:,:,bs)';
   
    for sc=1:3
        bscore = cogscores(:,sc);
        
        [BBcc,BBp] = corr(nFCc_vec_all,bscore,'rows','complete','type','s');
        
        temp = zeros(N); temp(inds) = squeeze(BBp); BBmat = temp+temp';
        mask = (BBmat<edge_thr).*~eye(N);
        [~,cpsz] = get_components(mask);
        cpsz_max = max(cpsz);
        mdx = module_density(mask,block_idx,1,'sum');
        mask_nnz = nnz(mask);
        
        % random permutations
        mdsum = zeros(nSys); cpszr_max = zeros(nPerm,1);
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
        
        % p-vals
        mds_pval = 1-(mdsum./nPerm);
        cpsz_pval = 1-(sum(cpsz_max>cpszr_max)/nPerm);
        
        % package
        
        ncc_out.modules_permP(:,:,bs,sc) = mds_pval;
        
        ncc_out.components_permP(bs,sc) = cpsz_pval;
        
        ncc_out.largest_components(bs,sc) = cpsz_max;
        
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
                ncc_out.corr_thr = mc_thr/((nSys^2-nSys)/2+nSys);
                ncc_out.plot(bs,sc) = nnz(triu(mds_pval<=ncc_out.corr_thr));
            case 'fdr'
                ncc_out.corr_type = multccorr;
                ncc_out.corr_thr(bs,sc) = FDR(mds_pval,mc_thr);
                ncc_out.plot(bs,sc) = nnz(triu(mds_pval<=ncc_out.corr_thr(bs,sc)));
            otherwise
                fprintf(2,'multccorr does not match existing options. Defaulting to Bonferoni.\n')
                ncc_out.corr_type = 'bonferroni';
                ncc_out.corr_thr = mc_thr/((nSys^2-nSys)/2+nSys);
                ncc_out.plot(bs,sc) = nnz(triu(mds_pval<=ncc_out.corr_thr));
        end
        
        if bs==1 && sc==1
            disp(['binset ','cscore ','nsig'])
        end
        disp(['  ',num2str(bs),'      ',num2str(sc),'     ',num2str(ncc_out.plot(bs,sc))])
    end
end