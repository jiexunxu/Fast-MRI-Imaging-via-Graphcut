function [V, VEEV, segIndices, V0coefs]=amalgamateSubspaceVectors(V0, VEEV0, segIndices, ...
    x0, sensitivities, conjSensitivities, G, data, objParams, nbrTermWeights, ...
    fullNbrMat, EsEdratio, LSIter, V0coefs, LSSmoothnessScaling, perturbScale, invokeLSModelPeriod, iter)    
    if mod(iter-1, invokeLSModelPeriod)==0        
        V0coefs=LSApproxModel(V0, x0, sensitivities, conjSensitivities, G, data, objParams, ...
            nbrTermWeights, fullNbrMat, EsEdratio, LSIter, LSSmoothnessScaling);
   %     V0coefs=zeros(size(V0, 2), 1);dicPerSeg=size(V0, 2)/length(segIndices);idx=mod(ceil(iter/30), dicPerSeg)+1;
   %     V0coefs(idx:dicPerSeg:length(V0coefs))=0.0666*0.1^ceil(iter/30);        
    end
    V0coefs=V0coefs.*(rand(length(V0coefs), 1)*perturbScale*2+1-perturbScale);

    disp('Start producing V and VEEV...');
    [V, VEEV]=generateSubspace(V0, VEEV0, V0coefs, segIndices);
end

function [V, VEEV]=generateSubspace(V0, VEEV0, V0coefs, segIndices)
    n=length(segIndices);
    S=zeros(size(VEEV0, 1), n);
    V=zeros(size(V0, 1), n);
    for i=1:n
        idx=segIndices{i}(1):segIndices{i}(2);  
        S(idx, i)=V0coefs(idx);
        V(:, i)=V0(:, idx)*V0coefs(idx);
    end
    S=sparse(S);
    VEEV=S'*VEEV0*S;
end
