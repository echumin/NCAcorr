function ncc_plot_edgecorr(ncc_out,cogLabels,sysLabels,bs,sc)

narginchk(3,5)

load('fcn/colormaps.mat','bluered_cmap')

[~,nSys,nBS,nSC] = size(ncc_out.modules_permP);

if ~exist('bs','var') || isempty(bs)
    bs = 1:1:nBS;
end
if ~exist('sc','var') || isempty(sc)
    sc = 1:1:nSC;
end
if length(cogLabels) ~= nSC
    fprintf(2,'Number of provided cogLabels does not match number of cog variables in ncc_out\n')
    return
end
if length(sysLabels) ~= nSys
    fprintf(2,'Number of provided system (Block) labels does not match number of blocks in ncc_out\n')
    return
end

for sc_idx = 1:length(sc)
    figure
    tiledlayout('flow')
    for bs_idx = 1:length(bs)
        switch ncc_out.corr_type
            case 'bonferroni'
                corr_thr = ncc_out.corr_thr;
            case 'fdr'
                corr_thr = ncc_out.corr_thr(bs(bs_idx),sc(sc_idx));
        end
        [rw,cl]=find(triu(ncc_out.modules_permP(:,:,bs(bs_idx),sc(sc_idx))<=corr_thr));
        [X, Y] = ncc_grid_blocks([rw,cl],ncc_out.block_idx);
        nexttile
        imagesc(ncc_out.BB_corr_mat(:,:,bs(bs_idx),sc(sc_idx))); axis square;
        hold on; plot(X,Y,'k','LineWidth',1.5)
        title(['Bin Set ' num2str(bs(bs_idx))]); 
        ncc_bLabels_onEdgemat(sysLabels,ncc_out.block_idx)
        colormap(bluered_cmap); colorbar
    end
    sgtitle({[cogLabels{sc(sc_idx)} ' Blocks p<=0.05 ' ncc_out.corr_type ' adjusted'],' '}); 
end