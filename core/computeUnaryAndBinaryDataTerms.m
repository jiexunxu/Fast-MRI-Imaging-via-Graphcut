function VtEtEV=computeUnaryAndBinaryDataTerms(V, sensit, G) 
%{ 
    Can't do because ND-sparse array is not supported
    V=reshape(V, [size(sensit, 1) size(sensit, 2) size(sensit, 3) size(V, 2)]);
    FSv=zeros([size(V, 1) size(V, 2) size(V, 3) size(sensit, 4) size(V, 4)]);
    for i=1:size(sensit, 4)
        FSv(:, :, :, i, :)=bsxfun(@times, V, sensit(:, :, :, i));
    end
    FSv=fft(FSv, [], 1);
    FSv=fft(FSv, [], 2);
    GFSv=bsxfun(@times, FSv, G);
    FSv=reshape(FSv, numel(sensit), size(V, 4));
    GFSv=reshape(GFSv, numel(sensit), size(V, 4));
    VtEtEV=FSv'*GFSv/(size(sensit, 1)*size(sensit, 2));
%}
   % numVoxels=size(sensit, 1)*size(sensit, 2)*size(sensit, 3);
   % V=reshape(V, numVoxels, numel(V)/numVoxels);
   
    VtEtEV=zeros(size(V, 2), size(V, 2));
    sizeX=[size(sensit, 1), size(sensit, 2), size(sensit, 3)];
    Vt=V';conjSensit=conj(sensit);numPixels=size(sensit, 1)*size(sensit, 2);
    parfor i=1:size(V, 2)
        ev=sensit;
        ev=bsxfun(@times, ev, reshape(full(V(:, i)), sizeX));
        ev=fft(ev, [], 1);
        ev=fft(ev, [], 2);
        ev=bsxfun(@times, ev, G);
        ev=fft(conj(ev), [], 1);
        ev=fft(ev, [], 2);
        ev=conj(ev)/numPixels;            
        ev=ev.*conjSensit;
    	eev=sum(ev, 4);
        VtEtEV(:, i)=Vt*double(eev(:));
    end
end
