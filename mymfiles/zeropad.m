function paddedB =  zeropad(B, aspectratio)
% B is a k-space image
% zero pads B so that its aspect ratio (y/x) is at most equal to aspectratio 

[m,n] = size(B);
if m/n < aspectratio 
    paddedB = B;    % do nothing
else
    nnew = round(m/aspectratio);
    del = round((nnew - n)/2);
    paddedB = [zeros(m, del), B, zeros(m, del)];
end;
    