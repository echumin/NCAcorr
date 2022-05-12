function ncc_bLabels_onEdgemat(sysLabels,block_idx) 

[~, If] = unique(block_idx,'first');
[~, Il] = unique(block_idx,'last');

for i=1:length(If)
    if i>1
        Imid(i) = (Il(i)-If(i))/2 + Il(i-1);
    elseif i==1
        Imid(i) = Il(i)/2;
    end
end

xticks(Imid); xticklabels(sysLabels)
yticks(Imid); yticklabels(sysLabels)