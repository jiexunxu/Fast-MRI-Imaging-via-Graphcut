% Computes y=H'Hx for multiple coils and multiple slices. Both x and y
% are in image space
% Note: y is in the IMAGE SPACE, not k-space
% The dimension of x can either be size(H'H, 2), or size(H'H, 2)*2. In
% former case, everything proceeds as usual. In latter case, x is assumed
% to be a real vector representation of the complex image slices. In this
% case, output y will also be a real vector representing the complex y.

% Input y must be y=sensit, i.e assign the sensit values to y
function y=HtHx(x, sensit, conjSensit, G)
    y=sensit;
    y=bsxfun(@times, y, x);
    y=fft(y, [], 1);
    y=fft(y, [], 2);
    y=bsxfun(@times, y, G);
    y=fft(conj(y), [], 1);
    y=fft(y, [], 2);
    y=conj(y)/(size(y, 1)*size(y, 2));  
    y=y.*conjSensit;
    y=sum(y, 4);
end