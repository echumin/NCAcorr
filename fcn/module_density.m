function [Mden, Mden_in, Mden_bw, Mden_in_lo, Mden_bw_hi] = module_density(M,ci,diag,operation)
%
% For a given input matrix M, compute average within and between block
%   weights for the input partition ci. Dimentions of M must be the same
%   size as ci, with M(i,:) nodes corresponding community/affiliation in ci(i).
%
% INPUTS:
%   M   - An NxN adjacency matrix.
%   ci  - An affiliation vector size (N,1), with each element corresponding
%         to nodes in M.
%  diag - 1 if main diagonal is to be discounted (assumes modules on main
%         diagonal are square)
% operation - How are the block data reduced: 'mean', 'median', or 'sum'.
%
% OUTPUTS:
%   Mden - A matrix of block values as indexed by unique values in ci
%          (input partition).
%   Mden_in - Average within block value
%   Mden_bw - Average between block value
%
%   Mden_in_lo - lowest within block value
%   Mden_bw_hi - highest between block value
%
% Version history:
%   original - Olaf Sporns, Indiana University, ?
%   updated documentation and optional arguments - 
%       Evgeny Chumin, Indiana University, 2021
mnum = max(ci);
Mden = zeros(mnum);

for i=1:mnum
    iind = logical(ci==i);
    for j=1:mnum
        jind = logical(ci==j);
        Mij = M(iind,jind);
        if ((i==j)&&(diag==0))
            Mij = M(iind,jind);
        end
        if ((i==j)&&(diag==1))
            Mij = M(iind,jind);
            mask = find(ones(size(Mij)).*~eye(size(Mij,1)));
            Mij = Mij(mask);
        end
        switch operation
            case 'mean'
                Mden(i,j) = nanmean(Mij(:));
            case 'median'
                Mden(i,j) = nanmedian(Mij(:));
            case 'sum'
                Mden(i,j) = nansum(Mij(:));
        end
    end
end

Mden(isnan(Mden)) = 0;

m1mask = eye(mnum); ff1 = find(m1mask);
m2mask = ones(mnum).*~eye(mnum); ff2 = find(m2mask);
% mean density of within/between module blocks
% note: modules composed of a single node may have density 0
Mden_in = nanmean(Mden(ff1));
Mden_bw = nanmean(Mden(ff2));
% lowest within vs highest between
Mden_in_lo = min(nonzeros(Mden(ff1)));
Mden_bw_hi = max(Mden(ff2));


