function out = mymedian(x, strel)
% strel is binary mask, assumed to be odd size
if nargin < 2,  strel = ones(3,3);  end
[nr, nc] = size(x);
off = round((size(strel)-1)/2);
cr = off(1);
cc = off(2);
tmp = x;

for i=1+cr:nr-cr
    for j=1+cc:nc-cc
        tt = x(i-cr:i+cr, j-cc:j+cc);
        q = tt(strel==1);
        tmp(i,j) = median(q);
    end
end
out = tmp;