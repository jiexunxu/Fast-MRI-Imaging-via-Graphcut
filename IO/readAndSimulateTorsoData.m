function [data, sensit]=readAndSimulateTorsoData(filename, G)  
    data=readTorsoData(filename);
    [data, sensit]=undersampleData(data, G);   
    data=single(data);sensit=single(sensit);
end

function data=readTorsoData(filename)
    data=read_Pfile_with_hdr(filename, 2, 'ieee-le');  
end


function [data, sensit]=undersampleData(data, G)
    [~, ~, nS, nCoils]=size(data);        
    sensit=zeros(size(data));   
    if isempty(G)
        G=ones(size(data, 1), size(data, 2));
    end
    for i=1:nS
        for j=1:nCoils
            slice=data(:, :, i, j);
            idx=find(abs(slice)==max(abs(slice(:))));
            [mx, my]=ind2sub(size(slice), idx);
            cx=round(size(slice, 1)/2);cy=round(size(slice, 2)/2);
            data(:, :, i, j)=circshift(slice, [cx-mx cy-my]);
        end
    end
    Gsensit=maskingMatrix('sensitivity', size(G, 1), size(G, 2));
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