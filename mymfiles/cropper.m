function croppedB =  cropper(B, sz)
% latest modif - id sz is 2vector, crops to size sz
% B is a k-space image
% crops rows of B so that its aspect ratio (y/x) is at most equal to aspectratio 
[m,n] = size(B);
if length(sz)==1
    aspectratio = sz;
	if m/n < aspectratio 
        croppedB = B;    % do nothing
	else
        mnew = round(n*aspectratio);
        del = round((m-mnew)/2);
        croppedB = B(del+1:end-del, :);
	end;
else
    mnew = sz(1);  nnew = sz(2);
    delr = round((m-mnew)/2);
    delc = round((n-nnew)/2);
    
	croppedB = B(delr+1:delr+mnew, delc+1:delc+nnew);

end
