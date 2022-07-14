function [FCcomponents, block_idx,dis] = fcn_binRSS(ts,sys_aff,binSets)
%%                              RSS Binning                              %%
% INPUTS
%   ts      - 3D Time by Node by Subject matrix
%   sys_aff - System affiliation vector (e.g. yeo7 labels)
%   binSets - 2D binsets by number of bins matrix of ones and zeros, where
%            each row is a set of bins: 1=included 0=excluded root sum squared
%            (RSS) ordered data bin.
%            [1 1 1 1] -> include all data (equivalent to pearson FC)
%            [1 0 0 0] -> include only time points from lowest 25% RSS 
%            (first of four bins). 
%
% OUPUTS
%   FCcomponents -> Data subset components based on binSets
%   block_idx    -> Node system affiliations that are consistent with
%                   reordered time series.
%   dis          -> nubmer of time points separation between adjacent
%                   points within RSS bins.

[block_idx,sys_ord] = sort(sys_aff,'ascend');

% Time by Node by Subject Time series
[T, ~, ~] = size(ts);

% size of bins
[nSets,nBins] = size(binSets);
binsz = floor(T/nBins);

% bin percentiles
tmp = 0:100/nBins:100;
l1 = tmp(1:end-1);
l2 = tmp(2:end);
clear tmp
for i=1:length(l1)
    blabels{i}=[num2str(l1(i)) '-' num2str(l2(i)) '%'];
end

%crop timeseries on both ends to match bin size*number of bins
t2c = (T-(binsz*nBins));
if rem(t2c,2)==0 % even
    t2c = t2c/2;
    ts = ts(t2c+1:end-t2c,:,:);
else
    t2c = floor(t2c/2);
    ts = ts(t2c+2:end-t2c,:,:);
end

[T, ~, S] = size(ts);

bID = zeros(T,1);
for b=1:nBins
    bID((b-1)*binsz+1:b*binsz) = b;
end

figure
y1 = linspace(.5,1.5,S)';

disp('Computing compoments for subject:')
for s=1:S    
    disp(num2str(s))
    zts = zscore(squeeze(ts(:,sys_ord,s)));
    ets = fcn_edgets(zts);
    rss = nansum(ets.^2,2).^0.5;            % nansum will ignore nan edges
    [~,rss2] = sort(rss,'ascend');          % smallest to largest
    
    for i=1:nBins
        tp = rss2(binsz*(i-1)+1:binsz*i,1);
        stp = sort(tp);
        ax = scatter(stp,zeros(binsz,1)+y1(s)+(i-1),'x');
        if s == 1
            cxmap(i,:)=ax.CData;
        else
            ax.CData = cxmap(i,:);
        end
        hold on
     %   plot([min(stp) max(stp)],zeros(2,1)+y1(s)+(i-1),'Color',ax.CData)
        dis(:,i,s)=diff(stp);
    end  
    
    for bs=1:nSets
        id = find(binSets(bs,:)==1);
        frames = rss2(logical(sum(bID==id,2)));
        FCcomponents(:,s,bs) = nanmean(ets(frames,:),1); %nanmean will ignore contribution of nan edges
    end
end
xlabel('Time Series')
ylabel('RSS Bins')
yticks(1:1:nBins)
yticklabels(blabels)
title('Subject time points that fell in a particular RSS bin')
