function imgout = resample_imgvol(img, oldvoxres, newvoxres, vol_class)
% resamples img of voxel size oldvoxres into new image of voxel size newvoxres
nblocks = 4;
[nrows, ncols, nslices] = size(img);
nfact = newvoxres./oldvoxres;

x = 1:nrows;
xi = 1:nfact(1):nrows;
y = 1:ncols;
yi = 1:nfact(2):ncols;
z = 1:nslices;
zi = 1:nfact(3):nslices;

if nslices>50
    imgout = zeros(length(xi), length(yi), length(zi), vol_class);
    delsl = round(nslices/nblocks);
    beginsl = 1; endsl = delsl; newsl = 0;
    lastblock = false;
    while ~lastblock
        img1 = img(:,:,beginsl:endsl);
        img1 = internal_resample(img1, oldvoxres, newvoxres, vol_class);
        imgout(:,:, newsl+1:newsl+size(img1,3)) = img1;
        newsl = newsl+size(img1,3);
        beginsl = beginsl+delsl;
        endsl = endsl+delsl;
        if endsl>nslices
            endsl = nslices;
            delsl = endsl-beginsl+1;
        end
        if beginsl>nslices
            lastblock = true;
        end
    end    
else
    imgout = internal_resample(img, oldvoxres, newvoxres, vol_class);
end    

if nargin>3 && ~isempty(vol_class)
   imgout = eval([vol_class '(imgout)']);
end


function img = internal_resample(img, oldvoxres, newvoxres, vol_class)
nfact = newvoxres./oldvoxres;

if ndims(img)==2
    [nr, nc] = size(img);
    x = 1:nr;
    xi = 1:nfact(1):nr;
    y = 1:nc;
    yi = 1:nfact(2):nc;
    q = double(img);
    q = interp1(x, q, xi, 'linear');
    q = q.';
    q = interp1(y, q, yi, 'linear');
    q = q.';
    if nargin>3 & ~isempty(vol_class)
       img = eval([vol_class '(q)']);
    else
       img = q;
    end
    return;
end

[nr, nc, ns] = size(img);

% first smooth volume per axial slice
% NEED TO FIX CONV OPERATION IN NON-DOUBLE
hblur = fspecial('gaussian',[5, 5], 2.3);
for sl = 1:ns
   q = img(:,:,sl);
   q = conv2(single(q), hblur, 'same');
   if nargin>3 && ~isempty(vol_class)
       img(:,:,sl) = eval([vol_class '(q)']);
   else
       img(:,:,sl) = q;
   end
end

x = 1:nr;
xi = 1:nfact(1):nr;
y = 1:nc;
yi = 1:nfact(2):nc;
z = 1:ns;
zi = 1:nfact(3):ns;
if isinteger(img)
    szimg = size(img);
    if length(xi) > length(x)
        img = cat(1, img, zeros(length(xi) - length(x), szimg(2), szimg(3) ));
    end
    szimg = size(img);
    if length(yi) > length(y)
        img = cat(2, img, zeros(szimg(1), length(yi) - length(y), szimg(3)));
    end
    szimg = size(img);
    if length(zi) > length(z)
        img = cat(3, img, zeros(szimg(1), szimg(2), length(zi) - length(z)));
    end
    [xx,yy] = meshgrid(x,y);
    [xxi,yyi] = meshgrid(xi,yi);
    for k = 1:ns
        q = img(1:nr, 1:nc, k);
        qq = interp2(xx, yy, single(q), xxi, yyi, 'linear');
        img(1:length(xi),1:length(yi),k) = eval([class(img) '(qq)']);
    end
    img = img(1:length(xi),1:length(yi), :);
    for i = 1:length(xi)
        for j = 1:length(yi)
            q = squeeze(img(i,j,1:ns));
            qq = interp1(z, single(q), zi, 'linear'); 
            img(i,j,1:length(zi)) = eval([class(img) '(qq)']);
        end
    end
    img = img(1:length(xi),1:length(yi), 1:length(zi));
else
    img = interp1(x, img, xi, 'linear');
    img = permute(img, [2, 3, 1]);
    img = interp1(y, img, yi, 'linear');
    img = permute(img, [2, 3, 1]);
    img = interp1(z, img, zi, 'linear');
    img = permute(img, [2, 3, 1]);
end

if nargin>3 && ~isempty(vol_class)
   img = eval([vol_class '(img)']);
end

end

end