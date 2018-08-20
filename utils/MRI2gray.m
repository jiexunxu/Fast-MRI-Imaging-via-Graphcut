function img=MRI2gray(img, voxelMask, slice)
    img=abs(img);
    maxClr=max(img(:));minClr=min(img(:));    
    img=img-minClr;
    img=img/maxClr;
    [nR, nC, nS]=size(voxelMask);
    img=reshape(img, nR, nC, nS);
    img=img(:, :, slice);
    img=img/max(abs(img(:)));
end