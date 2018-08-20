function newsensit=processSensitMap(sensit, voxelMask, d) 
    [nr, nc, ns, ncl]=size(sensit);
    theta=angle(sensit);
    newsensit=zeros(nr, nc, ns, ncl);
    for k=1:ncl
        for j=1:ns
            newsensit(:, :, j, k)=imdilate(abs(sensit(:, :, j, k)), strel('disk', d, 0));
            newsensit(:, :, j, k)=imfilter(abs(newsensit(:, :, j, k)), fspecial('gaussian', [d, d], 3.0));
            newsensit(:, :, j, k)=newsensit(:, :, j, k).*voxelMask(:, :, j);
        end
    end
    newsensit=newsensit.*exp(1i*theta);
%{
    [nr, nc, ns, ncl]=size(sensit);
    newsensit=zeros(nr, nc, ns, ncl);
    diskmask=diskMask(d);
    for k=1:ncl
        for j=1:ns
            sensit1=sensit(:, :, j, k);
            [r, c, ~]=find(voxelMask(:, :, j));
            for i=1:length(r)
                x=r(i);y=c(i);
                mask=iterMask(x, y, nr, nc, diskmask);
                locals=sensit1.*mask;
                [~, idx]=max(abs(locals(:)));
                newsensit(x, y, j, k)=locals(idx);
            end
        end
    end
%}
end

function mask=iterMask(x, y, nr, nc, diskmask)
    mask=zeros(nr, nc);
    d=(size(diskmask, 1)-1)/2;    
    mask(x-d:x+d, y-d:y+d)=diskmask;
end

function mask=diskMask(d)
    mask=zeros(2*d+1, 2*d+1);
    for i=1:size(mask, 1)
        for j=1:size(mask, 2)
            if sqrt((i-(d+1))^2+(j-(d+1))^2)<=d
                mask(i, j)=1;
            end
        end
    end
end