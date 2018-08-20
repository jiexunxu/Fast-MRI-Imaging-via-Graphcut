function [roimask] = my2droidraw_old(img, refine_flg)
% general purpose 2d mask selection, allows multiple rois to be selected
% supports color images

nreps = 6;
if nargin<2, refine_flg = 0; end
[nrows, ncols, ncolors] = size(img);
zoomfact = 600/max(nrows, ncols);
q = img/max(img(:))*255;
roimask = zeros(nrows, ncols);

slice_done = false;
figure;
disp('draw your ROI polygon - use as many points as needed, right click when done. Esc or Enter for next slice');
while ~slice_done
    imagesc(uint8(q)); truesize(round([nrows, ncols]*zoomfact)); h = gcf; if ncolors==1, colormap(gray); end
    [bw, xi, yi] = roipoly;
    figure(h);
    if isempty(xi), 
        %disp('is empty');
        slice_done = true;
    else
        roimask(bw) = 1;
    end
end   
if nnz(roimask)==0,
    return;
end
%figure; imagesc(roimask); colormap(gray); title('user selected roi')
close(h);
%% refine selected ROI to conform to loval edges
% if nargin >1 && refine_flg
if refine_flg
    sig_edge = 0.8;
    h = fspecial('sobel');
    hblur = fspecial('gaussian', [5 5], sig_edge);
    roiperim = bwperim(roimask);
    roidil = bwmorph(roimask, 'dilate', 2);
    roidist = bwdist(roiperim);
    %roier = bwmorph(roimask, 'erode', 5);
    roier = (roidist>5 & roimask);
    stats = regionprops(roimask, 'MinorAxisLength');
    sigdist = 0.1*stats.MinorAxisLength;
    roidist_thr = 1.5*sigdist;
    edgemap = zeros(nrows,ncols);
    edg = zeros(nrows,ncols);
    for cl = 1:ncolors
        q = img(:,:,cl);
        qq = q;
        %qq = conv2(q, hblur, 'same');
        edgh = imfilter(qq, h, 'same');
        edgv = imfilter(qq, h.', 'same');
        edgemap = edgemap + sqrt(edgh.*edgh + edgv.*edgv);
        edg = edg + edge(q, 'canny', [], sig_edge);    
    end
    edg(roidist > roidist_thr) = 0;
    edg(roidil==0) = 0;
    edg = single(edg>=ncolors/2);
    edgemap(edg==0) = 0;
    edgemap = edgemap/max(edgemap(:)) + exp(-roidist.^2 / 2/ sigdist^2);
    edgevec = edgemap(:);
    edgevec = edgevec(edgevec>0);
    edgemap_thr = 0.5*median(edgevec);
    edg(edgemap < edgemap_thr) = 0;

   %% connect end points of separate edgesegments
    edg = bwareaopen(edg, 5);    
    %edg = bwmorph(edg, 'dilate');
    edg(roidil==0) = 0;
    %edg = imfill(edg, 'holes');
    for rep = 1:nreps
         edg = connect_edges(edg);
    end
    edg = imfill(edg, 'holes');
    edg = bwmorph(edg, 'dilate', 2);
    edg = imfill(edg, 'holes');
    if nnz(edg) < 0.6*nnz(roimask)
        edg1 = bwmorph(edg, 'skel', 4);
        for rep = 1:nreps
             edg1 = connect_edges(edg1);
        end
        edg1 = imfill(edg1, 'holes');
        edg1 = bwmorph(edg1, 'dilate', 2);
        edg1 = imfill(edg1, 'holes');
        edg = edg | edg1;
    end
    edg = edg | roier;
    %edg = bwmorph(edg, 'erode', 2);
    edg(roimask==0) = 0;
    roimask = edg;
    figure; imagesc(roimask); colormap(gray); title('refined roi')  
end

%% internal function
function edg = connect_edges(edgin)
        edg = edgin;
        %edg = bwmorph(edg, 'skel', 3);
        %edg = bwmorph(edg, 'remove');
        [L, ncomps] = bwlabel(edg);
        segends = [];
        for i = 1:ncomps
            edgeseg = single(L==i);
            tmp = edgeseg - bwmorph(edgeseg, 'spur');
            [ii, jj] = find(tmp);  % assumes only 2 ends
%             segends = [segends; [ii, jj]];
            if length(ii)<2
                segends = [segends; [ii, jj; ii, jj]];
            else
                Q = zeros(length(ii));
                for k=1:length(ii)
                    Q(:,k) = (ii-ii(k)).^2 + (jj-jj(k)).^2;
                end
                [maxd, kk] = max(Q(:));
                [comp1, comp2] = ind2sub(size(Q), kk);        
                segends = [segends; [ii([comp1, comp2]), jj([comp1,comp2])]];
            end
        end
        na = size(segends,1);
        A1 = repmat(segends(:,1) , [1,na]) - repmat((segends(:,1)).', [na,1]);
        A2 = repmat(segends(:,2) , [1,na]) - repmat((segends(:,2)).', [na,1]);
        A = sqrt(A1.*A1 + A2.*A2);
%         for i = 1:na,  A(i,i) = Inf; end
        for i = 1:2:na-1,  A(i:i+1,i:i+1) = Inf; end
        for i = 1:na
            i1 = segends(i,1); j1 = segends(i,2);
            [mind, k] = min(A(i,:));
            i2 = segends(k,1); j2 = segends(k,2);
            grad = (j2-j1)/(i2-i1+1e-3);
            mag = sqrt((i2-i1)^2 + (j2-j1)^2);
            % need to replace below with a better line drawing...
            edg1 = zeros(size(edg));
            for rho = linspace(0,mag, 100)
                ii = round(i1+rho/mag*(i2-i1)); ii = min(ii, nrows); ii = max(1, ii);
                jj = round(j1+rho/mag*(j2-j1)); jj = min(jj, ncols); jj = max(1, jj);
                edg1(ii,jj) = 1;
            end                
            edg1 = bwmorph(edg1, 'close');
            edg(edg1) = 1;
        end
         %edg = imfill(edg, 'holes');
        %edg = bwmorph(edg, 'close',2);
         %edg = bwmorph(edg, 'skel');
end % of internal fn

end % of main fn