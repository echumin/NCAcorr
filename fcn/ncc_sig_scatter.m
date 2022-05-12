function ncc_sig_scatter(FCcomponents,sig_block_idx,ncc_out,cogLabels,sysLabels,group,sig_flag,sign_flag)

narginchk(5,8)

[~,nSys,nBS,nSC] = size(ncc_out.modules_permP);

if length(cogLabels) ~= nSC
    fprintf(2,'Number of provided cogLabels does not match number of cog variables in ncc_out\n')
    return
end
if length(sysLabels) ~= nSys
    fprintf(2,'Number of provided system (Block) labels does not match number of blocks in ncc_out\n')
    return
end
[sBS,sSC] = size(sig_block_idx);
if sBS ~= nBS
    fprintf(2,'sig_blocks and nnc_out number of binsets do not match\n')
    return
end
if sSC ~= nSC
    fprintf(2,'sig_blocks and ncc_out number of cogScores do not match\n')
    return
end
if ~exist('sig_flag','var')
    sig_flag = 0;
end
if ~exist('sign_flag','var')
    sign_flag = 0;
end

cmap(1,:)=[0 0.4470 0.7410];
cmap(2,:)=[0.8500 0.3250 0.0980];
cmap(3,:)=[0.9290 0.6940 0.1250];
cmap(4,:)=[0.4940 0.1840 0.5560];

scmap = cmap(group,:);

N = length(ncc_out.block_idx);
[u,v] = find(triu(ones(N),1));  % get edges
u1 = sub2ind([N N],u,v); clear u v
% find significant binsets (r1) and cogscore pairs (c1)
[r1, c1] = find(~cellfun(@isempty,sig_block_idx));

for i=1:length(r1)
    BBcm = ncc_out.BB_corr_mat(:,:,r1(i),c1(i)); % pull out corr matrix
    
    % isolate first level significant edges
    switch sig_flag
        case 1
            uncPmask = logical(ncc_out.edge_uncorrPthr_mask(:,:,r1(i),c1(i)));
            sigmat = BBcm./uncPmask; sigmat(isinf(sigmat))=0;
        case 0
            sigmat = BBcm;
    end
    % find significant blocks
    [r2,c2]=ind2sub([nSys nSys],sig_block_idx{r1(i),c1(i)});
    fcc = FCcomponents(:,:,r1(i)); 
    cog = ncc_out.cogscores(:,c1(i)); 
    for j=1:length(r2)
        sbm = (ncc_out.block_idx==r2(j))*(ncc_out.block_idx==c2(j))';
        switch sig_flag
            case 1
                sbm = triu(and(sbm,uncPmask),1);
            case 0
                sbm = logical(triu(sbm,1));
            otherwise
                fprintf(2,'Unknown sig_flag! Plotting average of all in block edges\n')
                sbm = logical(triu(sbm,1));
        end
        
        figure
        
        switch sig_flag
            case 1
                st = '[avg-sig edges]';
            case 0
                st = '[avg-all edges]';
        end  
        
        switch sign_flag
            case 1
                psbm = and(sigmat>0,sbm);
                nsbm = and(sigmat<0,sbm);
                pos_sbm_edges = psbm(u1);
                neg_sbm_edges = nsbm(u1);
                pos_fcc_sig =  mean(fcc(pos_sbm_edges,:),1)';
                neg_fcc_sig =  mean(fcc(neg_sbm_edges,:),1)';
                
                subplot(1,2,1)
                scatter(pos_fcc_sig,cog)
                ylabel(cogLabels{c1(i)}); 
                xlabel({'Mean block FCcomp per subject','insert:edge corr distribution',[num2str(nnz(psbm)) ' edges']});
                title('pos corr edges')
                axes('Position',[.36 .22 .1 .2])
                boxplot(sigmat(psbm)); xticks([])
                subplot(1,2,2)
                scatter(neg_fcc_sig,cog)
                ylabel(cogLabels{c1(i)});
                xlabel({'Mean block FCcomp per subject','insert:edge corr distribution',[num2str(nnz(nsbm)) ' edges']});
                title('neg cor edges')
                axes('Position',[.8 .22 .1 .2])
                boxplot(sigmat(nsbm)); xticks([])
                sgtitle(sprintf('Binset %d: %s - %s block & %s %s',r1(i),sysLabels{r2(j)},sysLabels{c2(j)},cogLabels{c1(i)},st))
            case 0
                sbm_edges = sbm(u1);
                fcc_sig =  mean(fcc(sbm_edges,:),1)';
                [r,p]=corr(fcc_sig,cog,'rows','complete','type','Spearman');
                scatter(fcc_sig,cog,50,scmap,'filled')
                title(sprintf('Binset %d: %s - %s block & %s %s',r1(i),sysLabels{r2(j)},sysLabels{c2(j)},cogLabels{c1(i)},st))
                ylabel(cogLabels{c1(i)});
                xlabel({'Mean block FCcomp per subject','insert:edge corr distribution',['Block Avg: R=' num2str(r) ' p=' num2str(p)]});
                axes('Position',[.8 .22 .1 .2])
                boxplot(sigmat(sbm)); xticks([])
            otherwise
                fprintf(2,'Unknown sign_flag! Plotting positively and negatively correlated edges together\n')
                sbm_edges = sbm(u1);
                fcc_sig =  mean(fcc(sbm_edges,:),1)';
                [r,p]=corr(fcc_sig,cog,'rows','complete','type','Spearman');
                scatter(fcc_sig,cog)
                title(sprintf('Binset %d: %s - %s block & %s %s',r1(i),sysLabels{r2(j)},sysLabels{c2(j)},cogLabels{c1(i)},st))
                ylabel(cogLabels{c1(i)});
                xlabel({'Mean block FCcomp per subject','insert:edge corr distribution',['Block Avg: R=' num2str(r) ' p=' num2str(p)]});
                axes('Position',[.8 .22 .1 .2])
                boxplot(sigmat(sbm)); xticks([])
        end
    end
end