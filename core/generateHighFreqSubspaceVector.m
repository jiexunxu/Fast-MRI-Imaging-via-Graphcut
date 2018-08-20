function V=generateHighFreqSubspaceVector(x0, voxelMask, segVoxelIndices, BFsimgas, BFsigmaWeight, use3DBF, G)       
    [numRows, numCols, numSlices]=size(x0);   
    Vlen=0;
    for w=1:length(BFsimgas)
        Vlen=Vlen+length(BFsimgas{w});
    end
    V=zeros(numRows, numCols, numSlices, Vlen+4);
    ctr=1;
    for w=1:length(BFsimgas)
        for ictr=1:length(BFsimgas{w})
            sig=BFsimgas{w}(ictr);
            if use3DBF
                tmp=bilateralFiltering3D(x0, w, sig, BFsigmaWeight*sig)-x0;
            else
                tmp=bilateralFilter(x0, 2, [sig BFsigmaWeight*sig])-x0;
            end
            V(:, :, :, ctr)=tmp;
            ctr=ctr+1;
        end
    end
   
  
    W=swt2(x0, 4, 'db3');
    W2=zeros(size(W));    
    for i=2:size(W, 4)
    %    W2.dec{i}=wthresh(W.dec{i}, 's', 0.7*median(abs(W.dec{i}(:))));
        W2(:, :, :, i)=W(:, :, :, i);
        idx=floor((i+1)/3);
        V(:, :, :, Vlen+idx)=-V(:, :, :, Vlen+idx)-cast(iswt2(W2, 'db3'), 'double');
      %  V(:, :, i-1)=-indwt2(W2);
        W2(:, :, :, i)=zeros(size(W2(:, :, :, i)));
    end 
   % xtmp=ifft2(fftshift(fftshift(fft2(x0)).*G));
   % V(:, :, size(V, 4))=xtmp-x0;
    
    
    V=reshape(V, numel(x0), numel(V)/numel(x0));  
    V=sparsifyDenseSubspace(V, segVoxelIndices, voxelMask);
 %   V=V.*(0.9+0.2*rand(size(V)));
end
