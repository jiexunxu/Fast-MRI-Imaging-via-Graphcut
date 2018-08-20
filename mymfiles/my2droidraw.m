function [roimask] = my2droidraw(img, refine_flg)
% general purpose 2d mask selection, allows multiple rois to be selected
% supports color images

if nargin<2, refine_flg = 1; end
[nrows, ncols, ncolors] = size(img);
zoomfact = 800/max(nrows, ncols);
q = img/max(img(:))*255;
roimask = zeros(nrows, ncols);

slice_done = false;
%figure;
disp('draw your ROI polygon - use as many points as needed, right click when done. Esc or Enter for next slice');
while ~slice_done
%     subplot(1,2,2);
%     imagesc(zeros(size(q))); truesize(round([nrows, ncols]*zoomfact)); colormap(gray);
%     subplot(1,2,1);
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
%close(h);
%% refine selected ROI to conform to local edges
% if nargin >1 && refine_flg


nreps = 10;        % # of region growing steps
stdfact = 2.0;     % defines range of "body" intensities - # of stdevs from mean body intensity
sig_edge = 1.8;
h = fspecial('sobel');    
hblur = fspecial('gaussian', [5 5], sig_edge);

if refine_flg
%     roiperim = bwperim(roimask);
%     roidil = bwmorph(roimask, 'dilate', 2);
%     roidist = bwdist(roiperim);
%     %roier = bwmorph(roimask, 'erode', 5);
%     roier = (roidist>5 & roimask);
%     stats = regionprops(roimask, 'MinorAxisLength');
%     sigdist = 0.1*stats.MinorAxisLength;
%     roidist_thr = 1.5*sigdist;

    edgemap = zeros(nrows,ncols);
    edg = zeros(nrows,ncols);
    for cl = 1:ncolors
        q = img(:,:,cl);
        qq = q;
        qq = conv2(q, hblur, 'same');
        edgh = imfilter(qq, h, 'same');
        edgv = imfilter(qq, h.', 'same');
        edgemap = edgemap + sqrt(edgh.*edgh + edgv.*edgv);
        edg = edg + edge(q, 'canny', [], sig_edge);    
%         edg = edg + edge(q, 'zerocross', [], h);    
    end
    edg = single(edg>=ncolors/2);
    %figure; imagesc(edg); colormap(gray);
    edg = bwmorph(edg, 'dilate', 1);
    edg(roimask==1) = 0;
    distmap = bwdist(edg);
    edgemap(edg==0) = 0;
    if ncolors>1,    img = sqrt(sum(img.^2, 3)); end
    roivec = img(roimask>0);
    thrlo = mean(roivec) - stdfact*std(roivec);
    thrhi = mean(roivec) + stdfact*std(roivec);

%     edgemap = edgemap/max(edgemap(:)) + exp(-roidist.^2 / 2/ sigdist^2);
%     edgevec = edgemap(:);
%     edgevec = edgevec(edgevec>0);
%     edgemap_thr = 0.5*median(edgevec);
%     edg(edgemap < edgemap_thr) = 0;

   %% Now grow region
    stop_growing = 0;
    for k = 1:nreps
        % first remove boundary points that extend beyond edges
        if mod(k+1,2)==0
            perim = bwperim(roimask);
            [ii,jj] = find(perim);
            for i = 1:length(ii)
                win = edgemap(ii(i)-1:ii(i)+1,jj(i)-1:jj(i)+1);
                ww = roimask(ii(i)-1:ii(i)+1,jj(i)-1:jj(i)+1);
                win = win.*ww;
                if max(max(win)) > win(2,2)  || edg(ii(i), jj(i))==1
                    roimask(ii(i),jj(i)) = 0;
                end
            end
        end
        roimaskold = roimask;
        roimask = bwmorph(roimask, 'dilate');
        roimask(distmap<1.2) = 0;
%         roimask(edg==1) = 0;
        if mod(k,2)==0
            delroi = roimask - roimaskold;
            if nnz(delroi) < 10
                stop_growing = 1;
            end
            roimask = bwmorph(roimask, 'hbreak');
            roimask = bwmorph(roimask, 'open', 2);
            roimask = bwareaopen(roimask, round(0.1*nnz(roimask)));
            roimask = imfill(roimask, 'holes');
        end
        if stop_growing,    break;   end
    end
    roimask = bwmorph(roimask, 'open');
    roimask(img<thrlo | img>thrhi) = 0;
    roimask = imfill(roimask, 'holes');
    %subplot(2,1,2); 
    imagesc(roimask); colormap(gray); title('refined roi') ; pause(1); 
end

end % of main fn