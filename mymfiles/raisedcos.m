function raised =  raisedcos(x, del, dimflg)
% returns a windowing by raised cosine
% handles 1d or 2d matrices - 
% in latter case does windowing in one or both dimensions
% dimflg: 'r' windows along rows only, 'c' along cols only, 'rc' for both
% del is the transition bracket on each side, in fraction of window length
% if different del's for rows and cols, pass a 2d vector

siz = size(x);
nrows = siz(1);
if length(siz)>1, ncols = siz(2);  end

if strcmp(dimflg, 'r')
    r = (rcos1d(ncols, del(1))).';
    rr = ones(nrows, 1)*r;
elseif strcmp(dimflg, 'c')
    r = rcos1d(nrows, del(1));
    rr = r*ones(1, ncols);
elseif strcmp(dimflg, 'rc')
    if length(size(del))>1, delr = del(1); delc = del(2);  
    else, delr = del; delc = del; end
    r = (rcos1d(ncols, delr)).';
    c = rcos1d(nrows, delc);
    rr = c*r;
end
raised = x.*rr;


function r = rcos1d(n, delfrac)
del = round(delfrac*n);
r = ones(n,1);
q = (cos((1:del).'/del*pi/2)).^2;
r = [flipud(q); r(del+1:end-del); q];
