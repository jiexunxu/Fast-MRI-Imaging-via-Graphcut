% Computes y=Hx for multiple coils and multiple slices. Both x and y
% are in image space
% Note: y is in the IMAGE SPACE, not k-space
% The dimension of x can either be size(H'H, 2), or size(H'H, 2)*2. In
% former case, everything proceeds as usual. In latter case, x is assumed
% to be a real vector representation of the complex image slices. In this
% case, output y will also be a real vector representing the complex y.
function y=Hx(x, sensitivities, G)    
    [numRows, numCols, numSlices, numCoils]=size(sensitivities);   
    x=reshape(x, [numRows numCols numSlices]);
    y=zeros(numRows, numCols, numSlices, numCoils);
    for slice=1:numSlices
        for i=1:numCoils
            signal=x(:, :, slice).*sensitivities(:, :, slice, i);
            signal=fft2(signal);            
            signal=signal.*G;
            signal=ifft2(signal);
            y(:, :, slice, i)=signal;
        end
    end
    y=reshape(y, numRows*numCols*numSlices*numCoils, 1);    
end