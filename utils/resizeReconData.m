function reconOut=resizeReconData(recon)
    if iscell(recon)
        N=length(recon);
        reconOut=cell(N, 1);
        for i=1:N
            tmp=zeros(512, 512, size(recon{i}, 3));
            for j=1:size(tmp, 3)
                tmp(:, :, j)=resizem(recon{i}(:, :, j), [512 512], 'bicubic');
            end
            reconOut{i}=tmp;
        end   
    else
        reconOut=zeros(512, 512, size(recon, 3));
        for i=1:size(recon, 3)
            reconOut(:, :, i)=resizem(recon(:, :, i), [512 512], 'bicubic');
        end
    end
end