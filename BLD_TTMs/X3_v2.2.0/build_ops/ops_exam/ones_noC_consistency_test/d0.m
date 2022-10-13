function [r1] = d0(v);
% make a diagonal matrix with v down the main diagonal
  m = length(v(:));
  r1 = spdiags([v(:)],[0],m,m);
