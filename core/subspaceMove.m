% algorithm: fmin(fmin), genetic algorithm (ga), simulated annealing (sa)
% vecCount: [1, size(V, 2)], specifies how many null vectors to randomly select 
% iterCount: total number of iterations to apply options.algorithm
% objParams: a 8-vector penalizer of the objective:
% lambda1*|E(x0+V*alpha)-y|+lambda2*SmoothnessTerm,
% objParams(1:2)=[lambda1 lambda2, objParams(3)= Es_threshold decay rate, objParams(4)=initial Es_threshold in
% the smoothness term, objParams(5)=power in the smoothness term,
% objParams(6)=minimum Es_threshold rate, objParams(7) is the neighbor count
% option for smoothnessTerm (either 6, 18 or 26), objParams(8) is bound for fminbnd, or iteration
% count for fminsearch and sa, or max runtime for ga

% opts=struct('inputMatData', struct('data', data, 'sensit', sensit,
% 'voxelMask', voxelMask, 'G', G), 'subspaceMove', struct('x0', x(:), 'subspace', V, 'iterCount', 9, 'objParams', [1e+8 10 -0.1 0.1 0.75 0.1 6 500], 'dataTermMatrix', VEEV));

function x=subspaceMove(opts)
    data=opts.inputMatData.data;
    sensitivities=opts.inputMatData.sensit;
    G=opts.inputMatData.G;    
    V=opts.subspaceMove.subspace; 
    x0=opts.subspaceMove.x0;
    nbrList=opts.subspaceMove.nbrList;
    nbrTermWeights=opts.subspaceMove.nbrTermWeights;
    iterCount=opts.subspaceMove.iterCount;
    dataTermMatrix=opts.subspaceMove.dataTermMatrix;
    objParams=opts.subspaceMove.objParams;
    fullNbrMat=opts.subspaceMove.fullNbrMat;
    numSegs=opts.subspaceMove.numSegs;    
    LSModelIter=opts.subspaceMove.LSModelIter;
    LSSmoothnessScaling=opts.subspaceMove.LSSmoothnessScaling;
    perturbScale=opts.subspaceMove.perturbScale;
    invokeLSModelPeriod=opts.subspaceMove.invokeLSModelPeriod;
    
    conjSensitivities=conj(sensitivities);
    dataTerm=@(x) norm(reshape(Hx(x, G, sensitivities)-data, numel(data), 1), 2)^2;
    Es_truncation_factor=objParams(3);
    objParams(3)=max(abs(x0(:)))*Es_truncation_factor;
    smoothnessTerm=@(x) smoothnessEnergy(x(:), objParams(3), objParams(4), nbrList, nbrTermWeights);    
    dataEnergy=objParams(1)*dataTerm(x0);smoothEnergy=smoothnessTerm(x0);
    EsEdratio=objParams(2);
    objParams(2)=dataEnergy/smoothEnergy*EsEdratio;    
    vecsPerSeg=size(V, 2)/numSegs;
    segIndices=cell(numSegs, 1);
    for i=1:numSegs
        segIndices{i}=[(i-1)*vecsPerSeg+1 i*vecsPerSeg];
    end     
    V0coefs=0;
    x=x0;
    % bk=E'(E(x0+sum_k-1 delta_xi)-y), iteratively updated
    bk=HtHx(x0, sensitivities, conjSensitivities, G)-Htx(data, conjSensitivities, G);
    bk=double(2*bk(:));
    constData=full(Hx(x0, G, sensitivities)-data);constData=constData(:);
    constData=real(objParams(1)*(constData'*constData));
    vxlSpts=[];    
    Es_old=objParams(1)*dataTerm(x0)+objParams(2)*smoothnessTerm(x0);
    for iter=1:iterCount            
      %  disp('amalgamate vectors...'); tic
        [W, newDataTermMatrix, segIndices, V0coefs]=amalgamateSubspaceVectors(V, dataTermMatrix, segIndices, ...
            x, sensitivities, conjSensitivities, G, data, objParams, nbrTermWeights, ...
            fullNbrMat, EsEdratio, LSModelIter, V0coefs, LSSmoothnessScaling, perturbScale, invokeLSModelPeriod, iter);                    
     %   time=toc;
     %   fprintf('amalgamation time=%d\n', time);        
     %   tic
        b=real(bk'*W);b=b';
        
        if isempty(vxlSpts)
            [vxlWeights, vxlSpts]=max(W, [], 2);
            vxlWeights=full(vxlWeights);vxlSpts(vxlWeights==0)=-1;
        end
     %   time=toc;
    %    fprintf('midterm time=%d\n', time);
     %   tic;
        alpha=QPBO_solve(W, x, newDataTermMatrix, b, nbrList, objParams, nbrTermWeights, vxlSpts, dataTerm, smoothnessTerm, constData);
     %   time=toc;
     %   fprintf('QPBO time=%d\n', time);        
        dx=W*alpha;dx=reshape(dx, size(x0));
        xnew=x+dx;Es=objParams(1)*dataTerm(xnew)+objParams(2)*smoothnessTerm(xnew);
        if Es_old>Es
            Es_old=Es;x=xnew;           
            bk=bk+2*reshape(HtHx(dx, sensitivities, conjSensitivities, G), numel(x0), 1);  
            fprintf('iter#=%d, committed moves=%d, total moves=%d, delta_x=%d, |x|=%d\n', iter, nnz(alpha), numel(alpha), norm(W*alpha), norm(x(:)));
        %    constData=full(Hx(x, G, sensitivities)-data);constData=constData(:);
        %    constData=real(objParams(1)*(constData'*constData));
        else
            fprintf('iter#=%d, %d moves are rejected |x|=%d\n', iter, nnz(alpha), norm(x(:)));
        end
        stats{iter}=uint8(255*abs(x)/max(abs(x(:))));
    end    
end
