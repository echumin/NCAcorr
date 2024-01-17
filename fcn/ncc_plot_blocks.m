function [sig_block_idx] = ncc_plot_blocks(ncc_out,cogLabels,sysLabels,bs,sc)
%%                   NCC Plot Blocks                   %%
% INPUTS
%   ncc_out      -> Data structure output by the fcn_netcontcorr code 
%   cogLabels    -> Cell variable where each element is a text label for
%                   the behavioral domains correlated with FCc in the main code.
%   sysLabels    -> Cell variable where elements are text labels for the
%                   blocks used in the main code.
%   bs           -> A vector of indices for which component blocks/bins you
%                   want the blocks plotted.
%   sc           -> A vector of indices for which behaviors/cognitive domains
%                   you want the blocks plotted.
%
% OUTPUTS
%   sig_block_idx -> a cell value where elements are behavior-FCcomponent
%                   pair indices that contain an indxed of significant network blocks.
%
% Evgeny Jenya Chumin, 2023, Indiana University
%%
narginchk(3,5)

load('fcn/colormaps.mat','red_cmap')

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
        nexttile
        imagesc(ncc_out.modules_permP(:,:,bs(bs_idx),sc(sc_idx))); axis square; caxis([0 .05])
        title(['Bin Set ' num2str(bs(bs_idx))]); 
        yticks(1:1:nSys); yticklabels(sysLabels); xticks(1:1:nSys); xticklabels(sysLabels)
        hold on;

        [sROW, sCOL]=find(triu(ncc_out.modules_permP(:,:,bs(bs_idx),sc(sc_idx))<=corr_thr));
        sig_block_idx{bs_idx, sc_idx} = sub2ind([nSys nSys],sROW,sCOL);
        clear sROW sCOL

        spy(ncc_out.modules_permP(:,:,bs(bs_idx),sc(sc_idx))<=corr_thr,'k*')

        colormap(flipud(red_cmap)); colorbar
        xlabel([]);
        ylim([0.5 (nSys+.5)]); xlim([0.5 (nSys+.5)])
    end
    sgtitle({[cogLabels{sc(sc_idx)} ' *p<=' num2str(corr_thr) ' ' ncc_out.corr_type ' adjusted'],' '}); 
end

