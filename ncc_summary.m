function ncc_summary(ncc_out,cogLabels,sysLabels,group)

%%
% 20k foot view, a binset by cog score matrix showing number of corrected
% significance.
addpath('NCAcorr/fcn')

figure;
tiledlayout(1,3,'TileSpacing','none')
ax(1)=nexttile;
    imagesc(ncc_out.binset)
    yticks(1:1:size(ncc_out.binset,1))
    xticks(1:1:size(ncc_out.binset,2))
    title('RSS Bin Sets')
ax(2)=nexttile(2, [1 2]);
    imagesc(ncc_out.plot)
    xticks(1:1:3); xticklabels(cogLabels);
    yticks([])
    colorbar
colormap (ax(1),[1 1 1 ; 0 0 0])
ax(2).XAxisLocation = 'top';

rg = (0:1/max(max(ncc_out.plot)):1)';
rcmap = [ones(max(max(ncc_out.plot))+1,1),flipud(rg),flipud(rg)];
colormap(ax(2),rcmap)  % replace this colormap with something stepwise

sgtitle('Number of MultComp Corrected Significant Blocks')
%%
% plot block by block matrices of p-values, marking corrected significance
sig_block_idx = ncc_plot_blocks(ncc_out,cogLabels,sysLabels) % replace w/ corrcoeff (maybe?)

% Plot edge-level adjacency matrices of correlations, highlighting the
% significant blocks
ncc_plot_edgecorr(ncc_out,cogLabels,sysLabels,[],[])

% function that uses tiledlayout flow to plot scatter plots for all significant
% blocks
ncc_sig_scatter(ncc_out.FCcomponents,sig_block_idx,ncc_out,cogLabels,sysLabels,group,0)

% plot fc components
ncc_plot_fccomp(ncc_out,sysLabels)