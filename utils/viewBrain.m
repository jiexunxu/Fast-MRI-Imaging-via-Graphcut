function [x, x2]=viewBrain(x, x2, varargin)
   % x=resizeBrain(x);x2=resizeBrain(x2);
    if nargin==2
        viewFullBrain(x, x2);
    else
        slice=varargin{1};
        if isempty(slice)
            return;
        end
        viewSlice(x, x2, slice);
    end
end

function viewFullBrain(x, x2)
    pause_time=0.16; 
    
    h=figure;set(h, 'Position', [0 0 1500 1000]);
    colormap(gray);
    for slice=1:size(x, 1)        
        subaxis(2, 3, 1, 'S', 0.01, 'M', 0);imagesc(squeeze(x(slice, :, :)));axis off;
        subaxis(2, 3, 4, 'S', 0.01, 'M', 0);imagesc(squeeze(x2(slice, :, :)));axis off;
        subaxis(2, 3, 2, 'S', 0.01, 'M', 0);imagesc(squeeze(x(:, slice, :)));axis off;
        subaxis(2, 3, 5, 'S', 0.01, 'M', 0);imagesc(squeeze(x2(:, slice, :)));axis off;
        subaxis(2, 3, 3, 'S', 0.01, 'M', 0);imagesc(squeeze(x(:, :, slice)));axis off;
        subaxis(2, 3, 6, 'S', 0.01, 'M', 0);imagesc(squeeze(x2(:, :, slice)));axis off;
        pause(pause_time);
    end    
end

function viewSlice(x, x2, slice)
    h=figure;set(h, 'Position', [0 0 1500 1000]);
    colormap(gray);
    subaxis(2, 3, 1, 'S', 0.01, 'M', 0);imagesc(squeeze(x(slice, :, :)));axis off;
    subaxis(2, 3, 4, 'S', 0.01, 'M', 0);imagesc(squeeze(x2(slice, :, :)));axis off;
    subaxis(2, 3, 2, 'S', 0.01, 'M', 0);imagesc(squeeze(x(:, slice, :)));axis off;
    subaxis(2, 3, 5, 'S', 0.01, 'M', 0);imagesc(squeeze(x2(:, slice, :)));axis off;
    subaxis(2, 3, 3, 'S', 0.01, 'M', 0);imagesc(squeeze(x(:, :, slice)));axis off;
    subaxis(2, 3, 6, 'S', 0.01, 'M', 0);imagesc(squeeze(x2(:, :, slice)));axis off;
end

function x=resizeBrain(x)
    x=abs(x);
    a=size(x);
    if sum(a~=[512 512 512])==0
        return
    end
    tmpx=zeros(256, 256, 256);
    sx=round((256-a(1))/2)+1;sy=round((256-a(2))/2)+1;sz=round((256-a(3))/2)+1;
    tmpx(sx:sx+a(1)-1, sy:sy+a(2)-1, sz:sz+a(3)-1)=x;
 %   x=interp3(abs(tmpx), 2, 'spline');
    x=tmpx;
    x=uint8(256*x/max(x(:)));
end
