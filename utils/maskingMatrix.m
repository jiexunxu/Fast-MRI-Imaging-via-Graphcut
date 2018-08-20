% This function returns the numRows-by-numCols masking matrix G based on
% particular sampling scheme used
function G=maskingMatrix(scheme, numRows, numCols)
    if strcmp(scheme, 'full')
        G=ones(numRows, numCols);
    elseif strcmp(scheme, 'random')
        % Random sampling with full sampling in the centeral area
        G=(rand(numRows, numCols)<0.07);
        rangeX=floor(numRows/2)-floor(numRows/8):floor(numRows/2)+floor(numRows/8);
        rangeY=floor(numCols/2)-floor(numCols/8):floor(numCols/2)+floor(numCols/8);    
        G(rangeX, rangeY)=ones(length(rangeX), length(rangeY));
    elseif strcmp(scheme, 'exponentialWeightedRandom')
        center=[ceil(numRows/2); ceil(numCols/2)];
        G=zeros(numRows, numCols);
        for i=1:numRows
            for j=1:numCols
                dist=norm([i-center(1) j-center(2)], 2);
                if rand<1/(1.011^dist)
                    G(i, j)=1;
                end                
            end
        end
    elseif strcmp(scheme, 'linearWeightedRandom')
        center=[ceil(numRows/2); ceil(numCols/2)];
        G=zeros(numRows, numCols);
        distMax=norm([center(1) center(2)], 2);
        probs=[1 1 0.15 0.005 0.005 0.005 0.001 0.001];
        N=length(probs);
        for i=1:numRows
            for j=1:numCols
                dist=norm([i-center(1) j-center(2)], 2);
                probsIdx=floor(0.95*N*dist/distMax)+1;
                if rand<probs(probsIdx)
                    G(i, j)=1;
                end                                    
            end
        end
    elseif strcmp(scheme, 'cartesian')
        % Cartesian sampling, skipping every other row
        G=zeros(numRows, numCols);
        G(1:5:numRows, :)=1;
    elseif strcmp(scheme, 'sensitivity')
        centerX=floor(numRows/2);centerY=floor(numCols/2);
        G=zeros(numRows, numCols);
        radX=ceil(numRows/12);radY=ceil(numCols/12);
        centerQuad=zeros(radX, radY);
        for i=1:radX
            for j=1:radY
                if (i*i)/(radX*radX)+(j*j)/(radY*radY)<=1
                    centerQuad(radX-i+1, j)=1; 
                end
            end
        end
        G(centerX-radX+1:centerX, centerY+1:centerY+radY)=centerQuad;
        G(centerX-radX+1:centerX, centerY-radY+1:centerY)=centerQuad(:, radY:-1:1);
        G(centerX+1:centerX+radX, centerY+1:centerY+radY)=centerQuad(radX:-1:1, :);
        G(centerX+1:centerX+radX, centerY-radY+1:centerY)=centerQuad(radX:-1:1, radY:-1:1);
    elseif strcmp(scheme, 'poisson')
        G=poissonSampling2D(numRows, numCols);
    elseif strcmp(scheme, 'poissonVariant')
        G=maskingMatrix('sensitivity', numRows, numCols);
        xC=ceil(size(G, 1)/2);yC=ceil(size(G, 2)/2);
        [existsX, existsY]=ind2sub(size(G), find(G));
        dists=sqrt((existsX-xC).^2+(existsY-yC).^2);
        [~, idx]=max(dists);
        px=existsX(idx);py=existsY(idx);
        list{1}=[px py];
        listPtr=0;
        iter=1;k=30;
        acc_factor=3.9;
        acc_exp=1;
        while listPtr<length(list)
            listPtr=listPtr+1;
            p=list{listPtr};
            r=(acc_factor*(sqrt(1/2)*sqrt(((p(1)-xC)/xC)^2/2+((p(2)-yC)/yC)^2/2)))^acc_exp;
            rads=r+rand(k, 1)*r;thetas=rand(k, 1)*2*pi;
            candP=zeros(k, 2);
            candP(:, 1)=round(p(1)+rads.*cos(thetas));
            candP(:, 2)=round(p(2)+rads.*sin(thetas));            
            candP(:, 1)=max(min(candP(:, 1), size(G, 1)), 1);
            candP(:, 2)=max(min(candP(:, 2), size(G, 2)), 1);
            for j=1:k
                r=(acc_factor*(sqrt(1/2)*sqrt(((candP(j, 1)-xC)/xC)^2+((candP(j, 2)-yC)/yC)^2)))^acc_exp;                
                if checkCollision(candP(j, 1), candP(j, 2), G, r)
                    G(candP(j, 1), candP(j, 2))=1;
                    list{length(list)+1}=candP(j, :);
                end
            end
            iter=iter+1;        
        end
    else
        fprintf('Unrecognized sampling method. Return full sampling matrix G\n');
        G=ones(numRows, numCols);
    end
end

function G=poissonSampling2D(nR, nC)
    G=maskingMatrix('sensitivity', nR, nC);
    minDists=[1.5 2 2.5];
    boundRFactors=[0.4 0.6 1];
    boundFactors=[1/6 1/4 2/5 1/2];
    xC=[ceil(nR/2) ceil(nC/2)];
    
    G2=zeros(max(nR, nC), max(nR, nC));
    if nR>nC
        G2(:, round((nR-nC)/2)+1:nC+round((nR-nC)/2))=G;
    else
        G2(round((nC-nR)/2)+1:nR+round((nC-nR)/2), :)=G;
    end
    for i=1:length(minDists)
        xMin=max(xC(1)-ceil(nR*boundFactors(i)), 1);
        xMax=min(xC(1)+ceil(nR*boundFactors(i)), nR);
        yMin=max(xC(2)-ceil(nC*boundFactors(i)), 1);
        yMax=min(xC(2)+ceil(nC*boundFactors(i)), nC);
        boundR=boundRFactors(i)*sqrt(2)*max(nR, nC)/2;
        G2=possionSampling2DHelper(G2, xMin, xMax, yMin, yMax, boundR, minDists(i));
    end
    if nR>nC
        G=G2(:, round((nR-nC)/2)+1:nC+round((nR-nC)/2));
    else
        G=G2(round((nC-nR)/2)+1:nR+round((nC-nR)/2), :);
    end
end

function G=possionSampling2DHelper(G, xMin, xMax, yMin, yMax, boundR, r)
    xC=ceil(size(G, 1)/2);yC=ceil(size(G, 2)/2);
    % Generate first seeding point  

    while true
   %     p=[ceil(rand()*(xMax-xMin)+xMin) ceil(rand()*(yMax-yMin)+yMin)];
        randR=rand()*boundR;randTheta=rand()*2*pi;
        px=round(xC+randR*cos(randTheta));
        py=round(xC+randR*sin(randTheta));
        px=max(min(px, size(G, 1)), 1);
        py=max(min(py, size(G, 2)), 1);
        if checkCollision(px, py, G, r) 
            G(px, py)=1;
            break;
        end
    end    

    %{
    [existsX, existsY]=ind2sub(size(G), find(G));
    dists=sqrt((existsX-xC).^2+(existsY-yC).^2);
    [~, idx]=max(dists);
    px=existsX(idx);py=existsY(idx);
    %}
    list{1}=[px py];
    listPtr=0;
    iter=1;k=30;
    while listPtr<length(list)
        listPtr=listPtr+1;
        p=list{listPtr};
        rads=r+rand(k, 1)*r;thetas=rand(k, 1)*2*pi;
        candP=zeros(k, 2);
        candP(:, 1)=round(p(1)+rads.*cos(thetas));
        candP(:, 2)=round(p(2)+rads.*sin(thetas));
        dists=sqrt((candP(:, 1)-xC).^2+(candP(:, 2)-yC).^2);
        thetas=atan2(candP(:, 1)-xC, candP(:, 2)-yC);       
        dists(dists>boundR)=boundR;
        candP(:, 1)=max(min(round(xC+dists.*cos(thetas)), size(G, 1)), 1);
        candP(:, 2)=max(min(round(yC+dists.*sin(thetas)), size(G, 2)), 1);        
    %    candP(:, 1)=max(min(candP(:, 1), xMax), xMin);
    %    candP(:, 2)=max(min(candP(:, 2), yMax), yMin);
        for j=1:k
            if checkCollision(candP(j, 1), candP(j, 2), G, r)
                G(candP(j, 1), candP(j, 2))=1;
                list{length(list)+1}=candP(j, :);
            end
        end
        iter=iter+1;        
    end
end

function collide=checkCollision(px, py, G, r)
    [existsX, existsY]=ind2sub(size(G), find(G));
    collide=min((existsX-px).^2+(existsY-py).^2)>r^2;
end