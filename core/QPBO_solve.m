function [alpha, unaryTerms, binaryTerms, u2, b2]=QPBO_solve(V, x0, VEEV, b, nbrList, objParams, nbrTermWeights, vxlSpts, dataTerm, smoothnessTerm, constData)     
    tic;  
    N=size(V, 2);
    unaryTerms=zeros(N, 2);
    unaryTerms(:, 2)=objParams(1)*real(diag(VEEV));
    [rows, cols, vals]=find(triu(VEEV, 1));
    binaryTerms=zeros(length(rows)+N*N, 6);
    binaryTermsCounter=length(rows)+1;
    binaryTerms(1:length(rows), 1)=rows;
    binaryTerms(1:length(rows), 2)=cols;
    binaryTerms(1:length(rows), 6)=objParams(1)*2*real(vals);            
    unaryTerms(:, 2)=unaryTerms(:, 2)+objParams(1)*b;
    
   % time=toc;
  %  fprintf('unary terms built in time %d\n', time);
    vxlWeights=full(sum(V, 2));
  %  fprintf('start binary terms\n');
  %  tic;
    [unaryTerms, binaryTerms, binaryTermsCounter]=graphcut_binaryterms_mex(size(unaryTerms, 1), unaryTerms(:), binaryTermsCounter, ...
    size(binaryTerms, 1), binaryTerms(:), size(nbrList, 1), nbrList(:)-1, real(x0), imag(x0), size(V, 1), vxlSpts, real(vxlWeights), ...
    imag(vxlWeights), objParams(3), objParams(4), objParams(2), nbrTermWeights);
  %  time=toc;
  %  fprintf('binary terms built in %d\n', time);
    unaryTerms=reshape(unaryTerms, N, 2);
    binaryTerms=reshape(binaryTerms, numel(binaryTerms)/6, 6);
  %  binaryTerms=binaryTerms(1:binaryTermsCounter-1, :);
    binaryTerms(binaryTerms(:, 1)==0, :)=[];
 %   b0=binaryTerms;u0=unaryTerms;
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
%    fprintf('start qpbo solve\n');
    alpha=vgg_qpbo(unaryTerms', uint32(binaryTerms(:, 1:2)'), binaryTerms(:, 3:6)');    
 %   fprintf('finish qpbo solve\n');
    alpha=double(alpha);alpha(alpha<0)=0;
    
 %   testUnaryAndBinaryEnergies(u0, b0, alpha, x0, V, objParams, dataTerm, smoothnessTerm, constData);
end
