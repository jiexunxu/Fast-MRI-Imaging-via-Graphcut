function voxelMask=generateVoxelMask(x)
    [nx, ny, nz]=size(x);    
    voxelMask=false(nx, ny, nz);
    for i=1:nz
        mask=x(:, :, i)>pi/3.6;
        mask=bwmorph(mask, 'dilate', 10);
        mask=bwmorph(mask, 'erode', 10);
        mask=not(mask);
        voxelMask(:, :, i)=mask;
    end
    voxelMask([1:11 nx-10:nx], :, :)=false;
    voxelMask(:, [1:11 ny-10:ny], :)=false;
  %  voxelMask(:, :, [1:11 nz-10:nz])=false;
    
    for i=1:nx
        voxelMask(i, :, :)=bwmorph(squeeze(voxelMask(i, :, :)), 'fill', 3);
        voxelMask(i, :, :)=bwmorph(squeeze(voxelMask(i, :, :)), 'dilate', 2);
        voxelMask(i, :, :)=maxComponent(squeeze(voxelMask(i, :, :)));
    end
    for i=1:ny
        voxelMask(:, i, :)=bwmorph(squeeze(voxelMask(:, i, :)), 'fill', 3);
        voxelMask(:, i, :)=bwmorph(squeeze(voxelMask(:, i, :)), 'dilate', 2);
        voxelMask(:, i, :)=maxComponent(squeeze(voxelMask(:, i, :)));
    end
    for i=1:nz
        voxelMask(:, :, i)=bwmorph(voxelMask(:, :, i), 'fill', 3);
        voxelMask(:, :, i)=bwmorph(voxelMask(:, :, i), 'dilate', 2);
        voxelMask(:, :, i)=maxComponent(squeeze(voxelMask(:, :, i)));
    end
    voxelMask=maxComponent(voxelMask);
    
    for i=1:nx
        voxelMask(i, :, :)=bwmorph(squeeze(voxelMask(i, :, :)), 'fill', 3);
        voxelMask(i, :, :)=bwmorph(squeeze(voxelMask(i, :, :)), 'majority', 3);
        voxelMask(i, :, :)=bwmorph(squeeze(voxelMask(i, :, :)), 'dilate', 5);   
        voxelMask(i, :, :)=bwmorph(squeeze(voxelMask(i, :, :)), 'erode', 5);     
    end
    for i=1:ny
        voxelMask(:, i, :)=bwmorph(squeeze(voxelMask(:, i, :)), 'fill', 3);
        voxelMask(:, i, :)=bwmorph(squeeze(voxelMask(:, i, :)), 'majority', 3);
        voxelMask(:, i, :)=bwmorph(squeeze(voxelMask(:, i, :)), 'dilate', 5);
        voxelMask(:, i, :)=bwmorph(squeeze(voxelMask(:, i, :)), 'erode', 5);  
    end
    for i=1:nz
        voxelMask(:, :, i)=bwmorph(voxelMask(:, :, i), 'fill', 3);
        voxelMask(:, :, i)=bwmorph(voxelMask(:, :, i), 'majority', 3);
        voxelMask(:, :, i)=bwmorph(voxelMask(:, :, i), 'dilate', 5);     
        voxelMask(:, :, i)=bwmorph(voxelMask(:, :, i), 'erode', 5);   
    end
end

function Y=maxComponent(X)
    Y=false(size(X));
    CC=bwconncomp(X);
    bestMask=[];maxMask=0;
    for j=1:length(CC.PixelIdxList)
        mask2=CC.PixelIdxList{j};
        if length(mask2)>maxMask
            maxMask=length(mask2);
            bestMask=mask2;
        end
    end
    Y(bestMask)=true;     
end