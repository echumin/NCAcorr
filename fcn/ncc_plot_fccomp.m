function ncc_plot_fccomp(ncc_out,sysLabels,bs)

narginchk(2,3)

load('fcn/colormaps.mat','bluered_cmap')

[~,~,nBS] = size(ncc_out.FCcomponents);

if ~exist('bs','var') || isempty(bs)
    bs = 1:1:nBS;
end

N = length(ncc_out.block_idx);
[u,v] = find(triu(ones(N),1));  % get edges
u1 = sub2ind([N N],u,v); clear u v

figure
tiledlayout('flow')
for bs_idx = 1:length(bs)
    fcc = nanmean(ncc_out.FCcomponents(:,:,bs(bs_idx)),2);
    fccmat = zeros(N,N);
    fccmat(u1) = fcc; fccmat = fccmat+fccmat';
    nexttile
    imagesc(fccmat); axis square;
    title(['Bin Set ' num2str(bs(bs_idx))]); 
    ncc_bLabels_onEdgemat(sysLabels,ncc_out.block_idx)
    colormap(bluered_cmap); colorbar
    caxis([-max(abs(fcc)) max(abs(fcc))])
end
sgtitle('FC Components')
