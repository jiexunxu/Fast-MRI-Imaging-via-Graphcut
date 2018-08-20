function mywrite_movie_vol(movie_name, grayvol, labelvol, alpha)

% display  volume as axial slices
% movie_name should include full path, else will be saved to current dir
% if 2 vols are passed, adds up before writing movie
% the first vol is treated as grayscale and the second as a label img (will be converted to rgb)

% relative wieghts of gray vs label (only if both vols passed)
if nargin<4
    alpha = [0.8, 0.2];
end
if nargin > 2 && ~isempty(grayvol)  && ~isempty(labelvol)
    mxvol1 = max(grayvol(:));
    nslices = 1;
    ncolors = 1;
    szv1 = size(grayvol);
    nrows = szv1(1); ncols = szv1(2);
    if ndims(grayvol)>2,
        nslices = szv1(3);
    end
    if ndims(grayvol)>3,
        ncolors = szv1(4);
    end
    zoomfact = 600/max(nrows, ncols);
    if ~isequal(size(labelvol), szv1(1:3))
            disp('error in writing volume movie: 2 volumes have different sizes!');
    end
    figure; %set(fig,'DoubleBuffer','on');
    for sl = 1:nslices 
        labelimg = single(label2rgb(labelvol(:,:,sl)));
        q = single(squeeze(grayvol(:,:,sl,:)));
        if ncolors == 1, 
            q = cat(3, q, q, q); 
        elseif ncolors==2,
            q = cat(3, q, zeros(nrows, ncols)); 
        elseif ncolors>3
            q = q(:,:,1:3);
        end
        q = q/mxvol1*255;
        labelimg = labelimg/max(labelimg(:))*255;
        qq = alpha(1)*single(q) + alpha(2)*single(labelimg);
        imagesc(uint8(qq)); truesize(round([nrows, ncols]*zoomfact));
        F(sl) = getframe;
    end
elseif nargin > 2 && isempty(grayvol)
    [nrows, ncols, nslices] = size(labelvol);
    zoomfact = 600/max(nrows, ncols);
    figure; %set(fig,'DoubleBuffer','on');
    for sl = 1:nslices 
        labelimg = single(label2rgb(labelvol(:,:,sl)));
        labelimg = labelimg/max(labelimg(:))*255;
        qq = single(labelimg);
        imagesc(uint8(qq)); truesize(round([nrows, ncols]*zoomfact));
        F(sl) = getframe;
    end
elseif nargin <= 2 || isempty(labelvol)
    mxvol1 = max(grayvol(:));
    nslices = 1;
    ncolors = 1;
    szv1 = size(grayvol);
    nrows = szv1(1); ncols = szv1(2);
    if ndims(grayvol)>2,
        nslices = szv1(3);
    end
    if ndims(grayvol)>3,
        ncolors = szv1(4);
    end
    zoomfact = 600/max(nrows, ncols);
    figure; %set(fig,'DoubleBuffer','on');
    for sl = 1:nslices 
        q = single(squeeze(grayvol(:,:,sl,:)));
        if ncolors == 1, 
            q = cat(3, q, q, q); 
        elseif ncolors==2,
            q = cat(3, q, zeros(nrows, ncols)); 
        elseif ncolors>3
            q = q(:,:,1:3);
        end
        q = q/mxvol1*255;
        imagesc(uint8(q)); truesize(round([nrows, ncols]*zoomfact)); 
        F(sl) = getframe; 
    end
end
movie2avi(F, movie_name, 'fps', 1, 'quality', 98, 'compression', 'Cinepak');
