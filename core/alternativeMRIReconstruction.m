function x=alternativeMRIReconstruction(data, sensit, G, x0, voxelMask)
    objParams=[1 0.8 1 0.8];
          
    [data, sensit, G, x0, voxelMask, xBkgrd]=preprocessInput(data, sensit, G, x0, voxelMask);
    xtmp=x0(:);xtmp(xtmp==0)=[];objParams(3)=objParams(3)*median(abs(xtmp));   
    [segVoxelIndices, segmentedNbrLists, nbrList]=preprocessSegmentation(abs(x0), voxelMask);
         
    V0=generateCandidateDictionary(x0, voxelMask, segVoxelIndices); 
    disp('Finish bilateral subspace');
    [VEEV0, ytEV0]=computeVEEV_ytEV(V0, sensit, G, data, x0);    
    disp('Finish VEEV and ytEV');
    x=generalizedJumpMoves(x0, V0, VEEV0, ytEV0, objParams, segmentedNbrLists, data, sensit, G, nbrList);
    x=x+xBkgrd;
end

function [data, sensit, G, x0, voxelMask, xBkgrd]=preprocessInput(data, sensit, G, x0, voxelMask)
    if isempty(voxelMask)
        voxelMask=manualImageSegmentation(abs(xSENSE));
    end
    x0=double(x0);
    xBkgrd=x0.*(~voxelMask);
    x0=x0.*voxelMask;
    for i=1:size(data, 4)
        data(:, :, :, i)=double(data(:, :, :, i)).*voxelMask;
        sensit(:, :, :, i)=double(sensit(:, :, :, i)).*voxelMask;
    end
end

function [segVoxelIndices, segmentedNbrLists, nbrList]=preprocessSegmentation(x0, voxelMask)
    nbrList=computeNbrList(voxelMask, true(size(x0)), [1 1]);
 %   [segVoxelIndices, segLabelMap]=segmentImageAndMat(abs(x0), voxelMask, 0.002); 
    [segVoxelIndices, segLabelMap]=segmentInput(x0, voxelMask); 
    segmentedNbrLists=produceSegmentedNbrLists(nbrList, segVoxelIndices, segLabelMap);
end

function [segVoxelIndices, segLabelMap]=segmentInput(x0, voxelMask)
    segLabelMap=zeros(size(x0));
    if length(size(x0))==2
        [I, J]=ind2sub(size(x0), find(voxelMask));
        minX=min(I);maxX=max(I);
        minY=min(J);maxY=max(J);
        xSegs=10;ySegs=10;
        xStep=ceil((maxX-minX)/xSegs);
        yStep=ceil((maxY-minY)/ySegs);
        maxLbl=0;
        for x=minX:xStep:maxX
            for y=minY:yStep:maxY
                BB=[x min(x+xStep-1, maxX) y min(y+yStep-1, maxY)];
                localMask=voxelMask(BB(1):BB(2), BB(3):BB(4));
                if nnz(localMask)>0
                    maxLbl=maxLbl+1;
                    segLabelMap(BB(1):BB(2), BB(3):BB(4))=maxLbl.*localMask;
                end
            end
        end
    elseif length(size(x0))==3
        [I, J, K]=ind2sub(size(x0), find(voxelMask));
        minX=min(I);maxX=max(I);
        minY=min(J);maxY=max(J);
        minZ=min(K);maxZ=max(K);
        xSegs=6;ySegs=6;zSegs=6;
        xStep=ceil((maxX-minX)/xSegs);
        yStep=ceil((maxY-minY)/ySegs);
        zStep=ceil((maxZ-minZ)/zSegs);
        maxLbl=0;
        for x=minX:xStep:maxX
            for y=minY:yStep:maxY
                for z=minZ:zStep:maxZ
                    BB=[x min(x+xStep-1, maxX) y min(y+yStep-1, maxY) z min(z+zStep-1, maxZ)];
                    localMask=voxelMask(BB(1):BB(2), BB(3):BB(4), BB(5):BB(6));
                    if nnz(localMask)>0
                        maxLbl=maxLbl+1;
                        segLabelMap(BB(1):BB(2), BB(3):BB(4))=maxLbl.*localMask;
                    end
                end
            end
        end
    end
    segVoxelIndices=cell(maxLbl, 1);
    for i=1:maxLbl        
        segVoxelIndices{i}=find(segLabelMap==i);
    end    
end

function V=generateCandidateDictionary(x0, voxelMask, segVoxelIndices)       
    params{1}=linspace(1, 2, 2);  
      
    paramSetLength=0;
    for i=length(params)
        paramSetLength=paramSetLength+length(params{i});
    end
    V=zeros(numel(x0), paramSetLength); 
    ctr=1;
    for w=1:length(params)
        for ictr=1:length(params{w})
            sigma=params{w}(ictr);
        %    v=bilateralFiltering3D(x0, w, sigma, 0.011*sigma);    
            v=bilateralFilter(x0, w, [sigma 0.011*sigma])-x0;     
            V(:, ctr)=v(:);
            ctr=ctr+1;            
        end
    end  
    V=sparsifyDenseSubspace(V, segVoxelIndices, voxelMask);
end

function [VEEV, ytEV]=computeVEEV_ytEV(V, sensit, G, data, x0)
    VEEV=zeros(size(V, 2), size(V, 2));
    ytEV=zeros(1, size(V, 2));
    sizeX=[size(sensit, 1), size(sensit, 2), size(sensit, 3)];
    Vt=V';conjSensit=conj(sensit);numPixels=size(sensit, 1)*size(sensit, 2);  
    yt=Hx(x0, G, sensit);yt=2*(yt(:)-data(:))';
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
        ytEV(i)=real(yt*ev(:));
    end
end

function xk=generalizedJumpMoves(x0, V0, VEEV0, ytEV0, objParams, segmentedNbrLists, data, sensit, G, nbrList)
    graphcutIterCount=20;  
    xk=x0(:);alpha_xk=zeros(size(V0, 2), 1);
    for i=1:graphcutIterCount
        alpha=0.1*ones(size(V0, 2), 1);
        [unaryTerms, binaryTerms]=constructUnaryAndBinaryTerms(xk, alpha_xk, alpha, V0, VEEV0, ytEV0, objParams, segmentedNbrLists, G, sensit, data, nbrList);
        beta=QPBOSolve(unaryTerms, binaryTerms);
        alpha_xk=alpha_xk+alpha.*beta;        
        xk=x0(:)+V0*alpha_xk;
        fprintf('commited moves=%d, total moves=%d, commit percentage=%d\n', nnz(beta), length(beta), nnz(beta)/length(beta));             
    end
    xk=reshape(xk, size(x0));
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

function [unaryTerms, binaryTerms]=constructUnaryAndBinaryTerms(xk, alpha_xk, alpha, V0, VEEV0, ytEV0, objParams, segmentedNbrLists, G, sensit, data, nbrList)
    numSegs=length(segmentedNbrLists.nbrListPerSeg);
    elementsPerSeg=size(V0, 2)/numSegs;   
    Es_threshold=objParams(3);
    power=objParams(4);
    Vbar=V0*diag(alpha);
    
    % Contributions from data fidelity term
    VEEV=diag(alpha)'*VEEV0*diag(alpha);
    ytEV=ytEV0*diag(alpha)+2*alpha_xk'*VEEV0*diag(alpha);
    constEd=Hx(reshape(xk, size(G, 1), size(G, 2)), G, sensit);constEd=constEd(:)-data(:);constEd=constEd'*constEd;
    unaryTerms=zeros(size(VEEV, 1), 2);        
    unaryTerms(:, 2)=unaryTerms(:, 2)+objParams(1)*real(diag(VEEV));
    unaryTerms(:, 2)=unaryTerms(:, 2)+objParams(1)*transpose(ytEV);
    [rows, cols, vals]=find(triu(VEEV, 1));
    binaryTerms=zeros(length(rows)+size(V0, 2)*elementsPerSeg*5, 6);
    binaryTerms(1:length(rows), 1)=rows;
    binaryTerms(1:length(rows), 2)=cols;
    binaryTerms(1:length(rows), 6)=objParams(1)*2*real(vals);
    binaryTermsCounter=length(rows)+1;
    
    uEd=unaryTerms;bEd=binaryTerms;
    % Contributions from smoothness term
    for i=1:numSegs
        segNbrList=segmentedNbrLists.nbrListPerSeg{i};
        segIdx=(i-1)*elementsPerSeg;
        for v1idx=1:elementsPerSeg
            unaryTerms(segIdx+v1idx, 1)=unaryTerms(segIdx+v1idx, 1)+objParams(2)*smoothnessEnergy(xk, Es_threshold, power, segNbrList);
            unaryTerms(segIdx+v1idx, 2)=unaryTerms(segIdx+v1idx, 2)+objParams(2)*smoothnessEnergy(xk+full(Vbar(:, segIdx+v1idx)), Es_threshold, power, segNbrList);
            for v2idx=v1idx+1:elementsPerSeg
                [binaryTerms, binaryTermsCounter]=updateBinaryTerms(binaryTerms, binaryTermsCounter, Vbar, segIdx+v1idx, segIdx+v2idx, xk, Es_threshold, power, segNbrList, objParams);
            end
        end
        for j=i+1:numSegs
            interNbrList=segmentedNbrLists.nbrListInterSeg{i, j};
            if isempty(interNbrList)
                continue;
            end
            seg2Idx=(j-1)*elementsPerSeg;
            for v1idx=1:elementsPerSeg
                for v2idx=1:elementsPerSeg
                    [binaryTerms, binaryTermsCounter]=updateBinaryTerms(binaryTerms, binaryTermsCounter, Vbar, segIdx+v1idx, seg2Idx+v2idx, xk, Es_threshold, power, interNbrList, objParams);
                end
            end            
        end
    end 
    uEs=unaryTerms-uEd;bEs=binaryTerms-bEd;
    dataTerm=@(x) norm(reshape(Hx(x, G, sensit)-data, numel(data), 1), 2)^2;
    smoothnessTerm=@(x) smoothnessEnergy(x(:), objParams(3), objParams(4), nbrList);    
    binaryTerms=binaryTerms(1:binaryTermsCounter-1, :);
    bEs=bEs(1:binaryTermsCounter-1, :);bEd=bEd(1:binaryTermsCounter-1, :);
end

function [binaryTerms, binaryTermsCounter]=updateBinaryTerms(binaryTerms, binaryTermsCounter, Vbar, v1idx, v2idx, xk, Es_threshold, power, segNbrList, objParams)
    v1=full(Vbar(:, v1idx));
    v2=full(Vbar(:, v2idx));
    binaryTerms(binaryTermsCounter, 1)=v1idx;
    binaryTerms(binaryTermsCounter, 2)=v2idx;
    % 00, 10, 01, 11
    binaryTerms(binaryTermsCounter, 3)=objParams(2)*smoothnessEnergy(xk, Es_threshold, power, segNbrList);
    binaryTerms(binaryTermsCounter, 4)=objParams(2)*smoothnessEnergy(xk+v1, Es_threshold, power, segNbrList);
    binaryTerms(binaryTermsCounter, 5)=objParams(2)*smoothnessEnergy(xk+v2, Es_threshold, power, segNbrList);
    binaryTerms(binaryTermsCounter, 6)=objParams(2)*smoothnessEnergy(xk+v1+v2, Es_threshold, power, segNbrList);
    binaryTermsCounter=binaryTermsCounter+1;
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

function E=smoothnessEnergy(x0, Es_threshold, power, nbrList)
    E=sum(min(abs(x0(nbrList(:, 1))-x0(nbrList(:, 2))).^power, Es_threshold));
end