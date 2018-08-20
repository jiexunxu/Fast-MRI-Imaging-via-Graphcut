function roimask = my3droidraw(vol)

refine_flg = 1;
nhood = ones(3,3,3); nhood = convn(nhood, nhood, 'same'); nhood = single(nhood>15);
strl = strel('arbitrary', nhood);

[nrows,ncols,nslices] = size(vol);
for sl = 1:nslices 
    disp('new slice');
    q = squeeze(vol(:,:,sl));
    roislice = my2droidraw(q, refine_flg);
    roimask(:,:,sl) = bwmorph(roislice, 'erode', 2);
    nnzvec(sl) = nnz(roislice);
    %fg = imfill(fg, 'holes');
end

% figure; title('ROI before 3d filtering');
% patch(isosurface(roimask,0.9), 'FaceColor',[0.75,.1,.75], 'EdgeColor','none');
% % view(45,30)
% %view(3)
% %axis tight
% daspect([1,1,1])
% lightangle(45,30);
% set(gcf,'Renderer','zbuffer'); lighting phong

ind = find(nnzvec>8);
% if min(ind)>1 
%     ind = [1:min(ind)-1, ind];
% end
% for sl = setdiff(1:nslices, ind)
for sl = setdiff(min(ind):max(ind), ind)
    %msk = zeros(rnows, ncols);
    i1 = max(ind(ind<sl));
    i2 = min(ind(ind>sl));
    alpha = (sl-i1) / (i2-sl);
    m1 = roimask(:,:, i1);
    m2 = roimask(:,:, i2);
    msk = single(m1 & m2);
    q = (m1 < m2);
    dist1 = bwdist(m1); 
    dist2 = bwdist(1-m2);   
    [ii,jj] = find(q);
    for i = 1:length(ii)
        if dist1(ii(i), jj(i)) <= alpha*dist2(ii(i),jj(i))
            msk(ii(i), jj(i)) = 1;
        end
    end
    q = (m1 > m2);
    dist1 = bwdist(1-m1); 
    dist2 = bwdist(m2);   
    [ii,jj] = find(q);
    for i = 1:length(ii)
        if dist1(ii(i), jj(i)) <= alpha*dist2(ii(i),jj(i))
            msk(ii(i), jj(i)) = 1;
        end
    end
    roimask(:,:, sl) = msk;
end

% finally, do 3d morph ops to clean up
roimask = imopen(roimask, strl);
roimask = imdilate(roimask, strl);

% display object in 3d
roimask = convn(double(roimask), ones(3,3,3), 'same');
roimask = roimask/max(roimask(:));

figure; title('User-drawn 3D ROI after 3d filtering');
patch(isosurface(roimask,0.8), 'FaceColor',[0.75,.1,.75], 'EdgeColor','none');
% view(45,30)
%view(3)
%axis tight
daspect([1,1,1])
lightangle(45,30);
set(gcf,'Renderer','zbuffer'); lighting phong

