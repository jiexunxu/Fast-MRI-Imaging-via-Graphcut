function visualize3DResult(varargin)
    N=length(varargin);
    h=figure;
    set(h, 'position', [0 0 1200 400*N]);
    colormap(gray);
    datas=varargin;
    visualize3DResultHelper(datas);
end

function visualize3DResultHelper(datas)
    [nx, ny, nz]=size(datas{1});
    N=max([nx ny nz]);
    M=length(datas);
    for j=1:N
        for i=1:M
            data=datas{i};
            x=min(j, nx);y=min(j, ny);z=min(j, nz);
            subaxis(M, 3, i*3-2, 'Margin', 0, 'Spacing', 0.02, 'Padding', 0);        
            imagesc(squeeze(abs(data(x, :, :))));
            subaxis(M, 3, i*3-1, 'Margin', 0, 'Spacing', 0.02, 'Padding', 0);        
            imagesc(squeeze(abs(data(:, y, :))));
            subaxis(M, 3, i*3, 'Margin', 0, 'Spacing', 0.02, 'Padding', 0);        
            imagesc(squeeze(abs(data(:, :, z))));
            pause(0.1);
        end
    end    
end