function [segVoxelIndices, L]=segmentImageAndMat(x0, voxelMask, ratio)

    [numRows, numCols, numSlices]=size(x0);    
    x=reshape(x0, numRows*numCols, numSlices);
    x=[x zeros(numRows*numCols, 1)];
    x=abs(x)/max(abs(x(:)));
    x=uint8(x*255);
    x=repmat(x, [1 1 3]);
    minSize=round(nnz(voxelMask)*ratio);
    L=vgg_segment_gb(x, 1.1, 1.1, minSize, numRows);
    L(:, size(L, 2))=[];
    L=reshape(L, size(x0));
    L(voxelMask==0)=0;
    allLabels=unique(L(:));allLabels(1)=[];    
    lbCount=length(allLabels);
    segVoxelIndices=cell(lbCount, 1);
    for i=1:lbCount
        L(L==allLabels(i))=i;
        segVoxelIndices{i}=find(L==i);
    end 
%{
    segLabelMap=zeros(size(x0));
    if length(size(x0))==2
        [I, J]=ind2sub(size(x0), find(voxelMask));
        minX=min(I);maxX=max(I);
        minY=min(J);maxY=max(J);
        xSegs=10;ySegs=10;
        xStep=ceil((maxX-minX)/xSegs);
        yStep=ceil((maxY-minY)/ySegs);
        maxLbl=0;
        for x=minX:xStep:maxX
            for y=minY:yStep:maxY
                BB=[x min(x+xStep-1, maxX) y min(y+yStep-1, maxY)];
                localMask=voxelMask(BB(1):BB(2), BB(3):BB(4));
                if nnz(localMask)>0
                    maxLbl=maxLbl+1;
                    segLabelMap(BB(1):BB(2), BB(3):BB(4))=maxLbl.*localMask;
                end
            end
        end
    elseif length(size(x0))==3
        [I, J, K]=ind2sub(size(x0), find(voxelMask));
        minX=min(I);maxX=max(I);
        minY=min(J);maxY=max(J);
        minZ=min(K);maxZ=max(K);
        xSegs=6;ySegs=6;zSegs=6;
        xStep=ceil((maxX-minX)/xSegs);
        yStep=ceil((maxY-minY)/ySegs);
        zStep=ceil((maxZ-minZ)/zSegs);
        maxLbl=0;
        for x=minX:xStep:maxX
            for y=minY:yStep:maxY
                for z=minZ:zStep:maxZ
                    BB=[x min(x+xStep-1, maxX) y min(y+yStep-1, maxY) z min(z+zStep-1, maxZ)];
                    localMask=voxelMask(BB(1):BB(2), BB(3):BB(4), BB(5):BB(6));
                    if nnz(localMask)>0
                        maxLbl=maxLbl+1;
                        segLabelMap(BB(1):BB(2), BB(3):BB(4))=maxLbl.*localMask;
                    end
                end
            end
        end
    end
    segVoxelIndices=cell(maxLbl, 1);
    for i=1:maxLbl        
        segVoxelIndices{i}=find(segLabelMap==i);
    end    
    %}
end

% Assume labels go from 1 to lbCount. labels(i, j, k)=0 indicate that the
% voxel is background
function [labels, lbCount]=segmentImage(x0, voxelMask, fraction)
    [numRows, numCols, numSlices]=size(voxelMask);
    labels=zeros(numRows, numCols, numSlices);
    lbCount=1;   
    for slice=1:numSlices
        img=MRI2gray(x0(:), voxelMask, slice);        
        img=round(img*255);
        tmp=zeros(numRows, numCols, 3);
        for j=1:3
            tmp(:, :, j)=img;
        end
        tmp=imfilter(tmp, fspecial('gaussian'));
        tmp=uint8(tmp);
        segSigma=1.1;
      %  segMinSize=round(nnz(voxelMask(:, :, slice))*fraction)+round((rand()-0.5)*nnz(voxelMask(:, :, slice))*fraction*0.1);    
        segMinSize=round(nnz(voxelMask(:, :, slice))*fraction);
        L=vgg_segment_gb(tmp, 0.8, segSigma, segMinSize);
        label=unique(L(:));
        for i=1:length(label)
            L(L==label(i))=lbCount;
            lbCount=lbCount+1;
        end        
        labels(:, :, slice)=L;
    end
    lbCount=lbCount-1;
    labels=labels.*voxelMask;
end