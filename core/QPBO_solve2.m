function [alpha, unaryTerms, binaryTerms, u2, b2]=QPBO_solve2(V, x0, VEEV, b, nbrList, objParams, nbrTermWeights, vxlSpts, dataTerm, smoothnessTerm, constData)     
    unaryTerms=zeros(size(V, 2), 2);
    unaryTerms(:, 2)=objParams(1)*real(diag(VEEV));
    [rows, cols, vals]=find(triu(VEEV, 1));
    binaryTerms=zeros(length(rows)+size(nbrList, 1), 6);
    binaryTermsCounter=length(rows)+1;
    binaryTerms(1:length(rows), 1)=rows;
    binaryTerms(1:length(rows), 2)=cols;
    binaryTerms(1:length(rows), 6)=objParams(1)*2*real(vals);        
    
    N=size(V, 2);
    unaryTerms(:, 2)=unaryTerms(:, 2)+objParams(1)*b;
    vxlWeights=full(sum(V, 2));
    [unaryTerms, binaryTerms, binaryTermsCounter]=graphcut_binaryterms_mex(size(unaryTerms, 1), unaryTerms(:), binaryTermsCounter, ...
    size(binaryTerms, 1), binaryTerms(:), size(nbrList, 1), nbrList(:)-1, real(x0), imag(x0), size(V, 1), vxlSpts, real(vxlWeights), ...
    imag(vxlWeights), objParams(3), objParams(4), objParams(2), nbrTermWeights);
    unaryTerms=reshape(unaryTerms, N, 2);
    binaryTerms=reshape(binaryTerms, numel(binaryTerms)/6, 6);
    binaryTerms(binaryTermsCounter:size(binaryTerms, 1), :)=[];
    
    % TODO: following code makes qpbo work, but makes energy not exact, need to get
    % most accurate way to do it
    unaryTerms(abs(unaryTerms)<1e-8*max(abs(unaryTerms(:))))=0;
    b2=binaryTerms(:, 3:6);
    b2(abs(b2)<1e-8*max(abs(b2(:))))=0;
    binaryTerms(:, 3:6)=b2;  
    u2=unaryTerms;
    u2(u2==0)=Inf;
    b2(b2==0)=Inf;
    minEnergy=(min([min(min(abs(u2))) min(min(abs(b2)))]));
    unaryTerms=unaryTerms/minEnergy*100000;
    binaryTerms(:, 3:6)=binaryTerms(:, 3:6)/minEnergy*100000;
    binaryTerms=int64(binaryTerms);
    unaryTerms=int64(unaryTerms);
    
    alpha=vgg_qpbo(unaryTerms', uint32(binaryTerms(:, 1:2)'), binaryTerms(:, 3:6)');    
    alpha=double(alpha);alpha(alpha<0)=0;
end
