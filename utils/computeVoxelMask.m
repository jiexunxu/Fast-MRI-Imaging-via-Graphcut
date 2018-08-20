function voxelMask=computeVoxelMask(x)
    weights=elipsoidWeights(size(x, 1), size(x, 2), size(x, 3));
    x=abs(angle(x));
    x=x.*weights;
    voxelMask=(x<pi/6);
    [m, n, p]=size(x);
    for i=1:p
        voxelMask(:, :, i)=bwmorph(voxelMask(:, :, i), 'clean');
        voxelMask(:, :, i)=bwmorph(voxelMask(:, :, i), 'spur');
        voxelMask(:, :, i)=bwmorph(voxelMask(:, :, i), 'close', 9);
        CC=bwconncomp(voxelMask(:, :, i));
        maxComponent=0;maxN=0;
        for j=1:CC.NumObjects
            N=length(CC.PixelIdxList{j});
            if N>maxN
                maxN=N;
                maxComponent=j;
            end
        end
        mask=false(m, n);
        mask(CC.PixelIdxList{maxComponent})=1;
        voxelMask(:, :, i)=mask;
        voxelMask(:, :, i)=bwmorph(voxelMask(:, :, i), 'close', 5);
        voxelMask(:, :, i)=bwmorph(voxelMask(:, :, i), 'diag', 5);
        voxelMask(:, :, i)=bwmorph(voxelMask(:, :, i), 'majority', 5);
    end
end

function weights=elipsoidWeights(numRows, numCols, numSlices)
    weights=zeros(numRows, numCols, numSlices);
    ratio10=4/5;ratio20=6/7;
    for i=1:numSlices
        shrinkFactor=1;ds=(abs(numSlices/2-i)/numSlices);
        ds=0.25;
        if ds>0.25
            shrinkFactor=1-(ds-0.25)*1.3;
        end
        ratio1=ratio10*shrinkFactor;ratio2=ratio20*shrinkFactor;
        a1=ratio1*numRows/2;b1=ratio1*numCols/2;
        a2=ratio2*numRows/2;b2=ratio2*numCols/2;
        weightsSlice=elipsoidWeightsSlice(numRows, numCols, a1, b1, a2, b2);
        weights(:, :, i)=weightsSlice;
    end
end

function weights=elipsoidWeightsSlice(numRows, numCols, a1, b1, a2, b2)
    center=[numRows/2 numCols/2];
    [rows, cols]=find(true(numRows, numCols));
    dx=rows-center(1);dy=cols-center(2);
    d1=(dx.^2)/(a1^2)+(dy.^2)/(b1^2);
    d2=(dx.^2)/(a2^2)+(dy.^2)/(b2^2);
    weights=ones(numRows, numCols);
    weights(d1>1 & d2<1)=1.2;
    weights(d2>1)=100;
end