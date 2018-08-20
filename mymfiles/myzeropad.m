function paddedB =  myzeropad(B, newsiz, unpad_fl, winflg)
% B is a k- or image-space image
% zero pads B symmetrically so that its new size is newsiz
% use unpad_fl=-1 to do unpaddin
% make sure the sizes are different by even number! - not any more :-)
% returns successfully even for odd differences, but result will probably
% be offset by 1 pixel
% Addition - now handles one-d vectors too
% Addition (4/13/05) - windowed if necessary - win is same size as paddedB

blur2d = [1 1 1; 1 1 1; 1 1 1];  blur2d = blur2d/sum(sum(blur2d));  % used for windowing
blur2d = conv2(blur2d, blur2d);
blur1d = [1 1 1].';  blur1d = blur1d/sum(sum(blur1d));
blur1d = conv2(blur1d, blur1d);

[m,n] = size(B);
if m==1      % B is a row vector
	mnew = 1;  nnew = newsiz;
    blur = blur1d.';
elseif n==1  % B is a col vector
    mnew = newsiz;  nnew = 1;
    blur = blur1d;
else         % B is a matrix
	mnew = newsiz(1);  nnew = newsiz(2);
    blur = blur2d;
end

if mnew==m & nnew==n 
    paddedB = B;    % do nothing
elseif nargin>2 & unpad_fl == -1    % un-zeropad to newsiz
    delr = round((m-mnew)/2+1)-1;
    delc = round((n-nnew)/2+1)-1;
    paddedB = B(delr+1:mnew+delr, delc+1:nnew+delc);
else
    delr = round((mnew - m)/2 + 1) - 1;
    delc = round((nnew - n)/2 + 1) - 1;
    %delc, delr, newsiz, size(B),
    paddedB = zeros(mnew, nnew);
    paddedB(delr+1:m+delr, delc+1:n+delc) = B;
end

% windowing is pretty bad right now....
if nargin>3 & winflg,  
    %paddedB =  raisedcos(paddedB, del, dimflg) % - figure out what to supply b4 use
    win = double(abs(paddedB) > 0);
    win = conv2(win, blur, 'same');
    paddedB = paddedB.*win;  
end
