function a = fcn_edgets(ts,zflag)
% Richard Betzel, Indiana University
if nargin == 1
    zflag = true;
end
[~,n] = size(ts);               % number samples/nodes
if zflag
    z = zscore(ts);                 % z-score
else
    z = ts;
end
[u,v] = find(triu(ones(n),1));  % get edges
a = z(:,u).*z(:,v);             % edge ts products