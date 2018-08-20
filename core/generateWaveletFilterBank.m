function V=generateWaveletFilterBank(x0, N, wname, voxelMask, segVoxelIndices) 
    V=nondecimatingWaveletMultiLevel(x0, N, wname);
    V=reshape(V, numel(x0), numel(V)/numel(x0));
    V=sparsifyDenseSubspace(V, segVoxelIndices, voxelMask);
end

function V=nondecimatingWaveletMultiLevel(x0, N, wname)
    if N==1
        thresholdFactors{1}=0.4:0.1:0.7;
    elseif N==2
        thresholdFactors{1}=0.0:0.1:0.4;
        thresholdFactors{2}=0:0.1:0.3;
    elseif N==3
        thresholdFactors{1}=0.2:0.1:0.5;
        thresholdFactors{2}=0.1:0.1:0.3;
        thresholdFactors{3}=0:0.1:0.2;
    elseif N==4
        thresholdFactors{1}=0.3:0.1:0.5;
        thresholdFactors{2}=0.1:0.1:0.3;
        thresholdFactors{3}=0:0.1:0.2;
        thresholdFactors{4}=0:0.1:0.2;
    end
    reconSize=0;
    for i=1:N
        reconSize=reconSize+length(thresholdFactors{i});
    end
    reconSize=reconSize*3;
    V=zeros(size(x0, 1), size(x0, 2), size(x0, 3), reconSize);
    
    %{
    if N>1
        error('Multi-level non-decimated 3D wavelet transform not yet implemented');
    end 
    coefsBase=nondecimatedWaveletTransform3D_mex(x0);
    N=size(coefsBase, 3)/8;
    for i=1:length(thresholdFactors{1})
        coefs=coefsBase;coefs(:, :, 1:N)=0;
        for j=1:7
            coef=coefs(:, :, j*N+1:(j+1)*N);
            thres=sqrt(2)*std(coef(:))*thresholdFactors{1}(i);
            coefH=my_wthresh(coef, 'h', thres);  
            coefs(:, :, j*N+1:(j+1)*N)=coefH;
        end
        v=nondecimatedWaveletReconstruction3D_mex(coefs);
        V(:, :, :, i)=-v;                    
    end
    %}
    for slice=1:size(x0, 3)        
        v=nondecimatingWaveletFilterBankRecursive2(squeeze(x0(:, :, slice)), wname, N, thresholdFactors, reconSize);
        V(:, :, slice, :)=-v;
    end
    
end

function recons=nondecimatingWaveletFilterBankRecursive2(img, wname, N, thresholdFactors, reconSize)
 %   W=ndwt2(img, N, wname);
    W=my_ndwt2(img, N);
    Wrecon=W;
    for i=1:length(Wrecon.dec)
        Wrecon.dec{i}=zeros(size(W.dec{i}));
    end
    reconsCtr=1;
    recons=zeros(size(img, 1), size(img, 2), reconSize);
    for i=2:length(Wrecon.dec)
        coefs=W.dec{i};        
        level=floor((length(Wrecon.dec)-i)/3)+1;
        for j=1:length(thresholdFactors{level})
            thres=sqrt(2)*std(coefs(:))*thresholdFactors{level}(j);
            coefsH=my_wthresh(coefs, 'h', thres);            
            Wrecon.dec{i}=coefsH;            
            recons(:, :, reconsCtr)=my_indwt2(Wrecon);
            reconsCtr=reconsCtr+1;
            Wrecon.dec{i}=zeros(size(W.dec{i}));
        end
    end
    recons(:, :, reconsCtr:size(recons, 3))=[];
end

function recons=nondecimatingWaveletFilterBankRecursive(img, wname, N, topLevel)
    if N<=0
        recons=[];
        return;
    end
    W=ndwt2(img, 1, wname);W2=W;reconsCtr=1;
    for i=1:4
        W2.dec{i}=zeros(size(W.dec{1}));
    end
    recons=zeros([size(img) power(3, N)+4]);
    if ~topLevel
        W2.dec{1}=W.dec{1};
        recons(:, :, reconsCtr)=indwt2(W2);reconsCtr=reconsCtr+1;        
        W2.dec{1}=zeros(size(W.dec{1}));
    end
    for i=2:length(W.dec)        
        W2.dec{i}=W.dec{i};
        recons(:, :, reconsCtr)=indwt2(W2);reconsCtr=reconsCtr+1;        
        W2.dec{i}=zeros(size(W.dec{1}));   
    end
    if ~topLevel
        reconsCA=nondecimatingWaveletFilterBankRecursive(W.dec{1}, wname, N-1, false);
        if ~isempty(reconsCA)
            for i=1:size(reconsCA, 3)
                W2.dec{1}=reconsCA(:, :, i);
                recons(:, :, reconsCtr)=indwt2(W2);reconsCtr=reconsCtr+1; 
                W2.dec{1}=zeros(size(W.dec{1})); 
            end
        end
    end
    reconsCH=nondecimatingWaveletFilterBankRecursive(W.dec{2}, wname, N-1, false);
    if ~isempty(reconsCH)
        for i=1:size(reconsCH, 3)
            W2.dec{2}=reconsCH(:, :, i);
            recons(:, :, reconsCtr)=indwt2(W2);reconsCtr=reconsCtr+1; 
            W2.dec{2}=zeros(size(W.dec{1})); 
        end
    end
    reconsCV=nondecimatingWaveletFilterBankRecursive(W.dec{3}, wname, N-1, false);
    if ~isempty(reconsCV)      
        for i=1:size(reconsCV, 3)
            W2.dec{3}=reconsCV(:, :, i);
            recons(:, :, reconsCtr)=indwt2(W2);reconsCtr=reconsCtr+1; 
            W2.dec{3}=zeros(size(W.dec{1})); 
        end
    end    
    recons(:, :, reconsCtr:size(recons, 3))=[];
end