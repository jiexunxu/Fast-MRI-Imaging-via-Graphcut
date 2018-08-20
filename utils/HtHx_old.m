% Computes y=H'Hx for multiple coils and multiple slices. Both x and y
% are in image space
% Note: y is in the IMAGE SPACE, not k-space
% The dimension of x can either be size(H'H, 2), or size(H'H, 2)*2. In
% former case, everything proceeds as usual. In latter case, x is assumed
% to be a real vector representation of the complex image slices. In this
% case, output y will also be a real vector representing the complex y.
function y=HtHx(x, sensitivities, G, voxelMaskList)    
    [numRows, numCols, numSlices, numCoils]=size(sensitivities);
    isreal=false;
    shrink_y=false;    
    if ~isempty(voxelMaskList) && length(x)==length(voxelMaskList)*2
        tmpx=zeros(numRows*numCols*numSlices*2, 1);
        tmpx([voxelMaskList*2-1; voxelMaskList*2])=x([1:2:length(x)-1 2:2:length(x)]);
        x=tmpx;
        shrink_y=true;
    elseif ~isempty(voxelMaskList) && length(x)==length(voxelMaskList)
        tmpx=zeros(numRows*numCols*numSlices, 1);
        tmpx(voxelMaskList)=x;
        x=tmpx; 
        shrink_y=true;
    elseif ~isempty(voxelMaskList) && length(x)==numRows*numCols*numSlices
        tmpx=zeros(numRows*numCols*numSlices, 1);
        tmpx(voxelMaskList)=x(voxelMaskList);
        x=tmpx;
    elseif ~isempty(voxelMaskList) && length(x)==numRows*numCols*numSlices*2 
        tmpx=zeros(numRows*numCols*numSlices*2, 1);
        tmpx([voxelMaskList*2-1; voxelMaskList*2])=x([voxelMaskList*2-1; voxelMaskList*2]);
        x=tmpx;
    elseif ~isempty(voxelMaskList)
        error('ApproximateNullVec:Hx', 'input x dimension is incorrect for the given voxelMaskList');
    end
    
    if length(x)==numRows*numCols*numSlices*2   
        isreal=true;
    elseif length(x)~=numRows*numCols*numSlices
        error('ApproximateNullVec:Hx', 'input x dimension is incorrect');
    end
    if isreal
        x=x(1:2:length(x)-1)+1i*x(2:2:length(x));
    end
    if issparse(x)
        x=full(x);
    end
    x=reshape(x, [numRows numCols numSlices]);
    y=zeros(numRows, numCols, numSlices);
    for slice=1:numSlices
        for i=1:numCoils
            signal=x(:, :, slice).*sensitivities(:, :, slice, i);
            signal=fft2(signal);
            signal=signal.*G;
            signal=ifft2(signal);
            y(:, :, slice)=y(:, :, slice)+signal.*conj(sensitivities(:, :, slice, i));
        end
    end
    y=reshape(y, numRows*numCols*numSlices, 1);
    
    if shrink_y
        y=y(voxelMaskList);
    end
    
    if isreal
        yreal=real(y);
        yimag=imag(y);
        y=zeros(length(y)*2, 1);
        y(1:2:length(y)-1)=yreal;
        y(2:2:length(y))=yimag;
    end
end