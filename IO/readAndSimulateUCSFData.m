function [data, sensit]=readAndSimulateUCSFData(filename, G)  
    if ischar(filename) && isempty(G)
        data=readUCSFData(filename);
        data=parseUCSFData(data, G); 
        return;    
    elseif ischar(filename)
        data=readUCSFData(filename);
        data=parseUCSFData(data, G); 
    else
        data=filename;
    end
    Gsensit=maskingMatrix('sensitivity', size(G, 1), size(G, 2));
    [data, sensit]=undersampleData(data, G, Gsensit);   
    data=single(data);sensit=single(sensit);
end

function data=readUCSFData(filename)
    data=readMultiRAID(filename);
    data=data{3}.image_obj.data;    
end

function [data, G]=parseUCSFData(data, G)
    data=data(1:2:size(data, 1), :, :, :);
    data=permute(data, [3 4 1 2]);
    
    [~, ~, nSlices, nCoils]=size(data);  
    % Zero pad data to 256-by-256
    %{
    if numRows>256 || numCols>256
        disp('Warning: Data has more than 256 rows or columns, no zero padding is done');
    else
        rowStartIdx=floor((256-numRows)/2)+1;
        colStartIdx=floor((256-numCols)/2)+1;
        data2=zeros(256, 256, nSlices, nCoils);
        data2(rowStartIdx:rowStartIdx+numRows-1, colStartIdx:colStartIdx+numCols-1, :, :)=data;
        data=data2;
        G2=zeros(256, 256);
        G2(rowStartIdx:rowStartIdx+numRows-1, colStartIdx:colStartIdx+numCols-1)=G;
        G=G2;
        clear data2 G2;
    end
    %}
    for i=1:nCoils
        data(:, :, :, i)=ifftshift(ifft(data(:, :, :, i), [], 3));        
    end    
    data=data(:, :, :, :);nSlices=size(data, 3);
    for i=1:nSlices
        for j=1:nCoils
            data(:, :, i, j)=fftshift(fft2(ifftshift(ifft2(data(:, :, i, j)))));
        end
    end
end

function [data, sensit]=undersampleData(data, G, Gsensit)
    [~, ~, nS, nCoils]=size(data);        
    sensit=zeros(size(data));   
    for i=1:nS
        for j=1:nCoils
            sensit(:, :, i, j)=ifft2(data(:, :, i, j).*G.*Gsensit);
            data(:, :, i, j)=ifft2(data(:, :, i, j).*G);
        end
    end
    sensit=processSensitMap(sensit, ones(size(data, 1) , size(data, 2), size(data, 3)), 13);    
    sensit=sensit/max(max(max(max(abs(sensit)))));
    data=data/max(max(max(max(abs(data)))));
end