function [x, xIndices]=fusionMoveMRIReconstruction(data, sensit, G, x0, voxelMask)
    objParams=[1 0.1 0.1 0.7];
    
    objParams(3)=objParams(3)*median(abs(x0(:)));  
    for i=1:size(data, 4)
        data(:, :, :, i)=data(:, :, :, i).*voxelMask;
        sensit(:, :, :, i)=sensit(:, :, :, i).*voxelMask;
    end
    [segVoxelIndices, segLabelMap]=segmentImageAndMat(abs(x0), voxelMask, 0.002);     
    nbrList=computeNbrList(voxelMask, true(size(x0)), [1 1]);
    
    V0=generateCandidateDictionary(x0, voxelMask, segVoxelIndices); 
    disp('Finish bilateral subspace');
    [VEEV0, ytEV0]=computeVEEV_ytEV(V0, sensit, G, data);    
    disp('Finish VEEV and ytEV');
    segmentedNbrLists=produceSegmentedNbrLists(nbrList, segVoxelIndices, segLabelMap);
    disp('Finish segmentedNbrLists');
    [x, xIndices]=fusionMoves(V0, VEEV0, ytEV0, objParams, segmentedNbrLists, data, sensit, G, nbrList);
    x=x+x0.*(~voxelMask);
end

function V=generateCandidateDictionary(x0, voxelMask, segVoxelIndices)       
    params{1}=linspace(2, 4, 10);  
      
    paramSetLength=0;
    for i=length(params)
        paramSetLength=paramSetLength+length(params{i});
    end
    V=zeros(numel(x0), paramSetLength+1); 
    V(:, 1)=x0(:);
    ctr=2;
    for w=1:length(params)
        for ictr=1:length(params{w})
            sigma=params{w}(ictr);
        %    v=bilateralFiltering3D(x0, w, sigma, 0.011*sigma);    
            v=bilateralFilter(x0, w, [sigma 0.011*sigma]);     
            V(:, ctr)=v(:);
            ctr=ctr+1;            
        end
    end  
    V=sparsifyDenseSubspace(V, segVoxelIndices, voxelMask);
end

function [VEEV, ytEV]=computeVEEV_ytEV(V, sensit, G, data)
    VEEV=zeros(size(V, 2), size(V, 2));
    ytEV=zeros(size(V, 2), 1);
    sizeX=[size(sensit, 1), size(sensit, 2), size(sensit, 3)];
    Vt=V';conjSensit=conj(sensit);numPixels=size(sensit, 1)*size(sensit, 2);    
    for i=1:size(V, 2)
        ev=sensit;
        ev=bsxfun(@times, ev, reshape(full(V(:, i)), sizeX));
        ev=fft(ev, [], 1);
        ev=fft(ev, [], 2);
        ev=bsxfun(@times, ev, G);
        ev=fft(conj(ev), [], 1);
        ev=fft(ev, [], 2);
        ev=conj(ev)/numPixels; 
        eev=ev.*conjSensit;        
    	eev=sum(eev, 4);
        VEEV(:, i)=Vt*double(eev(:));
        ytEV(i)=real(2*data(:)'*ev(:));
    end
end

function [x, xIndices]=fusionMoves(V0, VEEV0, ytEV0, objParams, segmentedNbrLists, data, sensit, G, nbrList)
    numSegs=length(segmentedNbrLists.nbrListPerSeg);
    candsPerSeg=size(V0, 2)/numSegs;
    fusionMovesIterCount=candsPerSeg*3;
    xIndices=ones(numSegs, 1);    
    for i=1:fusionMovesIterCount
        candIndices=ones(numSegs, 1)*(mod(i, candsPerSeg)+1);
        [unaryTerms, binaryTerms]=constructUnaryAndBinaryTerms(xIndices, candIndices, V0, VEEV0, ytEV0, objParams, segmentedNbrLists);
        alpha=QPBOSolve(unaryTerms, binaryTerms);
      %  nnz(alpha)
      %  dataTerm=@(x) objParams(1)*(norm(reshape(Hx(x, G, sensit)-data, numel(data), 1), 2)^2-abs(data(:)'*data(:)));
      %  smoothnessTerm=@(x) objParams(2)*smoothnessEnergy(x(:), objParams(3), objParams(4), nbrList);
     %   test(unaryTerms, binaryTerms, dataTerm, smoothnessTerm, V0, xIndices, candIndices);
        fprintf('commited moves=%d, total moves=%d, commit percentage=%d\n', nnz(alpha), length(alpha), nnz(alpha)/length(alpha));
        xIndices(alpha)=candIndices(alpha);        
    end
    xSubset=indices2subset(xIndices, candsPerSeg, numSegs);
    x=subset2x(V0, xSubset);
    x=reshape(x, size(data, 1), size(data, 2), size(data, 3));
end

function alpha=QPBOSolve(unaryTerms, binaryTerms)
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
    alpha(alpha<0)=0;alpha=logical(alpha);
end

function [unaryTerms, binaryTerms]=constructUnaryAndBinaryTerms(xIndices, candIndices, V0, VEEV0, ytEV0, objParams, segmentedNbrLists)
    N=length(xIndices);
    candsPerSeg=size(V0, 2)/N;
    xSubset=indices2subset(xIndices, candsPerSeg, N);
    candSubset=indices2subset(candIndices, candsPerSeg, N);    
    xVEEVx=VEEV0(xSubset, xSubset);
    cVEEVc=VEEV0(candSubset, candSubset);
    cVEEVx=VEEV0(candSubset, xSubset);
    x=subset2x(V0, xSubset);
    xCand=subset2x(V0, candSubset);
    Es_threshold=objParams(3);
    power=objParams(4);
    
    % Contributions from data fidelity term
    unaryTerms=zeros(N, 2);    
    unaryTerms(:, 1)=objParams(1)*(real(diag(xVEEVx))-ytEV0(xSubset));
    unaryTerms(:, 2)=objParams(1)*(real(diag(cVEEVc))-ytEV0(candSubset));    
    [rows, cols, vals]=find(triu(xVEEVx, 1));
    binaryTerms=zeros(length(rows)+round(N*N/2), 6);
    binaryTerms(1:length(rows), 1)=rows;
    binaryTerms(1:length(rows), 2)=cols;
    % 0, 0
    binaryTerms(1:length(rows), 3)=objParams(1)*2*real(vals);
    % 1, 0, values are in upper triangle of cVEEVx
    [~, ~, vals]=find(triu(cVEEVx, 1));    
    binaryTerms(1:length(rows), 4)=objParams(1)*2*real(vals);
    % 0, 1, values are in lower triangle of cVEEVx, or upper triangle of cVEEVx'
    [~, ~, vals]=find(triu(cVEEVx', 1)); 
    binaryTerms(1:length(rows), 5)=objParams(1)*2*real(vals);
    [~, ~, vals]=find(triu(cVEEVc, 1));   
    % 1, 1
    binaryTerms(1:length(rows), 6)=objParams(1)*2*real(vals);
    binaryTermsCounter=length(rows)+1;  
    % Contributions from smoothness term
    for i=1:N
        segNbrList=segmentedNbrLists.nbrListPerSeg{i};
        E=objParams(2)*smoothnessEnergy(x, Es_threshold, power, segNbrList);
        unaryTerms(i, 1)=unaryTerms(i, 1)+E;
        E=objParams(2)*smoothnessEnergy(xCand, Es_threshold, power, segNbrList);
        unaryTerms(i, 2)=unaryTerms(i, 2)+E;
        for j=i+1:N
            interNbrList=segmentedNbrLists.nbrListInterSeg{i, j};
            if isempty(interNbrList)
                continue;
            end
            binaryTerms(binaryTermsCounter, 1)=i;
            binaryTerms(binaryTermsCounter, 2)=j;          
            % seg1=0, seg2=0
            E=objParams(2)*smoothnessEnergy(x, Es_threshold, power, interNbrList);
            binaryTerms(binaryTermsCounter, 3)=E;
            % seg1=1, seg2=0
            interIndices=xIndices;interIndices(i)=candIndices(i);
            interSubset=indices2subset(interIndices, candsPerSeg, N);
            xInter=subset2x(V0, interSubset);
            E=objParams(2)*smoothnessEnergy(xInter, Es_threshold, power, interNbrList);
            binaryTerms(binaryTermsCounter, 4)=E;
            % seg1=0, seg2=1
            interIndices=xIndices;interIndices(j)=candIndices(j);
            interSubset=indices2subset(interIndices, candsPerSeg, N);
            xInter=subset2x(V0, interSubset);
            E=objParams(2)*smoothnessEnergy(xInter, Es_threshold, power, interNbrList);
            binaryTerms(binaryTermsCounter, 5)=E;
            % seg1=1, seg2=1
            E=objParams(2)*smoothnessEnergy(xCand, Es_threshold, power, interNbrList);
            binaryTerms(binaryTermsCounter, 6)=E;            
            binaryTermsCounter=binaryTermsCounter+1;
        end
    end 
    binaryTerms=binaryTerms(1:binaryTermsCounter-1, :);
end

function segmentedNbrLists=produceSegmentedNbrLists(nbrList, segVoxelIndices, segLabelMap)
    N=length(segVoxelIndices);
    nbrListPerSeg=cell(N, 1);
    nbrListPerSeg_Counter=ones(N, 1);
    nbrListInterSeg=cell(N, N);
    nbrListInterSeg_Counter=ones(N, N);
    for i=1:N
        nbrListPerSeg{i}=zeros(length(segVoxelIndices{i})*3, 2);
        for j=i+1:N
            nbrListInterSeg{i, j}=zeros(length(segVoxelIndices{i}), 2);
        end
    end
    
    for i=1:size(nbrList, 1)
        n1=nbrList(i, 1);n2=nbrList(i, 2);
        seg1=segLabelMap(n1);seg2=segLabelMap(n2);
        if seg1==seg2
            nbrListPerSeg{seg1}(nbrListPerSeg_Counter(seg1), 1)=n1;
            nbrListPerSeg{seg1}(nbrListPerSeg_Counter(seg1), 2)=n2;
            nbrListPerSeg_Counter(seg1)=nbrListPerSeg_Counter(seg1)+1;
        else
            if seg1>=seg2
                seg1=segLabelMap(n2);
                seg2=segLabelMap(n1);
            end
            nbrListInterSeg{seg1, seg2}(nbrListInterSeg_Counter(seg1, seg2), 1)=n1;
            nbrListInterSeg{seg1, seg2}(nbrListInterSeg_Counter(seg1, seg2), 2)=n2; 
            nbrListInterSeg_Counter(seg1, seg2)=nbrListInterSeg_Counter(seg1, seg2)+1;
        end
    end
    for i=1:N
        nbrListPerSeg{i}=nbrListPerSeg{i}(1:nbrListPerSeg_Counter(i)-1, :);
        for j=i+1:N
            nbrListInterSeg{i, j}=nbrListInterSeg{i, j}(1:nbrListInterSeg_Counter(i, j)-1, :);
        end
    end
    segmentedNbrLists.nbrListPerSeg=nbrListPerSeg;
    segmentedNbrLists.nbrListInterSeg=nbrListInterSeg;
end

function subset=indices2subset(indices, candsPerSeg, numSegs)
    subset=(1:candsPerSeg:candsPerSeg*numSegs)';
    subset=subset+indices-1;
end

function x=subset2x(V0, subset)
    x=full(sum(V0(:, subset), 2));
end






