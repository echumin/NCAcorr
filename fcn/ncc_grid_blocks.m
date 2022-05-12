function [X, Y] = ncc_grid_blocks(bPS,block_idx) 

[~, If] = unique(block_idx,'first');
[~, Il] = unique(block_idx,'last');

X = []; Y = [];
for i=1:size(bPS,1)
    x1 = [If(bPS(i,2))-.5 Il(bPS(i,2))+.5 Il(bPS(i,2))+.5 If(bPS(i,2))-.5 If(bPS(i,2))-.5 NaN]; 
    y1 = [If(bPS(i,1))-.5 If(bPS(i,1))-.5 Il(bPS(i,1))+.5 Il(bPS(i,1))+.5 If(bPS(i,1))-.5 NaN];
    X = [X, x1];
    Y = [Y, y1];
end
  