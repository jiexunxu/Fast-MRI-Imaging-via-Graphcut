% Computes y=H'*x
function y=Htx(x, sensitivities, G, varargin)
    [numRows, numCols, numSlices, numCoils]=size(sensitivities); 
    if nargin==4
        voxelMaskList=varargin{1};
    else
        voxelMaskList=[];
    end
    if issparse(x)
        x=full(x);
    end
    if length(size(x))>1
        x=reshape(x, numRows*numCols*numSlices*numCoils, 1);
    end
    isreal=false;
    if issparse(x)
        x=full(x);
    end
    if length(x)==numRows*numCols*numSlices*numCoils*2   
        isreal=true;
    elseif length(x)~=numRows*numCols*numSlices*numCoils
        error('ApproximateNullVec:Hx', 'input x dimension is incorrect');
    end
    if isreal
        x=x(1:2:length(x)-1)+1i*x(2:2:length(x));
    end
    x=reshape(x, numRows, numCols, numSlices, numCoils);
    y=zeros(numRows, numCols, numSlices);
    parfor slice=1:numSlices
        for i=1:numCoils            
            signal=fft2(x(:, :, slice, i));            
            signal=signal.*G;
            signal=ifft2(signal);
            y(:, :, slice)=y(:, :, slice)+signal.*conj(sensitivities(:, :, slice, i));            
        end
    end    
    y=reshape(y, numRows*numCols*numSlices, 1);
    if isreal
        yreal=real(y);
        yimag=imag(y);
        y=zeros(length(y)*2, 1);
        y(1:2:length(y)-1)=yreal;
        y(2:2:length(y))=yimag;
    end
    if ~isempty(voxelMaskList)        
        y=y(voxelMaskList);
    end
end