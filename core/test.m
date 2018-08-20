function test(uEd, bEd, uEs, bEs, Vbar, xk, dataTerm, smoothnessTerm, constEd)  
 %   beta=rand(size(Vbar, 2), 1)<0.5;
    beta=zeros(size(Vbar, 2), 1);
    Ed=getEnergies(beta, uEd, bEd)+constEd;
    Es=getEnergies(beta, uEs, bEs);
    Ed0=dataTerm(reshape(xk+Vbar*beta, 256, 176));
    Es0=smoothnessTerm(xk+Vbar*beta);
    fprintf('Ed=%d, Ed0=%d, EdErr=%d\n', Ed, Ed0, abs(Ed-Ed0)/Ed0);
    fprintf('Es=%d, Es0=%d, EsErr=%d\n', Es, Es0, abs(Es-Es0)/Es0);
end

function E=getEnergies(alpha, unaryTerms, binaryTerms)
    E=0;
    for i=1:size(unaryTerms, 1)
        E=E+unaryTerms(i, alpha(i)+1);                   
    end
    
    for i=1:size(binaryTerms, 1)
        if binaryTerms(i, 1)==0 || binaryTerms(i, 2)==0
            continue;
        end
        a1=alpha(binaryTerms(i, 1));a2=alpha(binaryTerms(i, 2));        
        if a1==0 && a2==0
            eng=binaryTerms(i, 3);            
        elseif a1==1 && a2==0
            eng=binaryTerms(i, 4);
        elseif a1==0 && a2==1
            eng=binaryTerms(i, 5);
        else
            eng=binaryTerms(i, 6);
        end
        E=E+eng;        
    end
end

%{
function test(v1, v2)
    fields=fieldnames(v1);
    for i=1:length(fields)
        mem1=v1.(fields{i});
        mem2=v2.(fields{i});
        if isstruct(mem1)
            test(mem1, mem2);
        elseif isnumeric(mem1)
            err=max(abs(mem1(:)-mem2(:)))/max(abs(mem1(:)));
            if err>1e-7
                disp(fields{i});
            end
        end
    end
end
%}
%{


function solIndices=mixIndices(alpha, xIndices, candIndices)
    solIndices=xIndices;
    solIndices(alpha==1)=candIndices(alpha==1);
end

function x=subset2x(V0, subset)
    x=full(sum(V0(:, subset), 2));
    x=reshape(x, 256, 176);
end

function subset=indices2subset(indices, candsPerSeg, numSegs)
    subset=(1:candsPerSeg:candsPerSeg*numSegs)';
    subset=subset+indices-1;
end

function [Ed, Es]=getEnergies(alpha, unaryTerms, binaryTerms)
    Ed=0;Es=0;N=length(alpha);
    for i=1:size(unaryTerms, 1)
        Ed=Ed+unaryTerms(i, alpha(i)+1);                   
    end
    for i=1:(N*N-N)/2
        a1=alpha(binaryTerms(i, 1));a2=alpha(binaryTerms(i, 2));        
        if a1==0 && a2==0
            eng=binaryTerms(i, 3);            
        elseif a1==1 && a2==0
            eng=binaryTerms(i, 4);
        elseif a1==0 && a2==1
            eng=binaryTerms(i, 5);
        else
            eng=binaryTerms(i, 6);
        end
        Ed=Ed+eng;        
    end
    
    for i=(N*N-N)/2+1:size(binaryTerms, 1)
        a1=alpha(binaryTerms(i, 1));a2=alpha(binaryTerms(i, 2));        
        if a1==0 && a2==0
            eng=binaryTerms(i, 3);            
        elseif a1==1 && a2==0
            eng=binaryTerms(i, 4);
        elseif a1==0 && a2==1
            eng=binaryTerms(i, 5);
        else
            eng=binaryTerms(i, 6);
        end
        Es=Es+eng;        
    end
end
%}

%{
function mask=test(figName, mask)   
    load(figName, '-mat');    
    figs=zeros(512, 512, 4);
    deltas{1}=0.58;alphas{1}=0.18;
    deltas{2}=0.6;alphas{2}=0.18;
    deltas{3}=0.59;alphas{3}=0.05;
    deltas{4}=1;alphas{4}=0.1;    
       
    for i=[4 1 2 3] 
       data=hgS_070000.children(i).children(1).properties.CData;              
       if size(data, 1)~=512
           data=resizem(data, [512 512], 'bicubic');
       end       
       if i==3
        %   data=imfilter(data, fspecial('gaussian', 2, 0.8)); 
       end
     %  data=data(512:-1:1, :);
       data=abs(data)/max(abs(data(:)));
       if isempty(mask)
           mask=manualImageSegmentation(abs(data));
       end
       if i==2
       %    data=imfilter(data, fspecial('gaussian', 2, 0.4));
       elseif i==3
      %     data=imfilter(data, fspecial('gaussian', 2, 0.8));
       end
       dataBrain=data.*mask;
       dataSkull=data.*(~mask);
       tmp=abs(dataBrain(:));tmp(tmp<1e-10)=[];dataBrainQ=max(tmp);
       tmp=abs(dataBrain(:));tmp(tmp<1e-10)=[];dataSkullQ=max(tmp);
       dataBrain=dataBrain/dataBrainQ;
       dataSkull=dataSkull/dataSkullQ;
   %    dataBrain=adapthisteq(dataBrain, 'NumTiles', [32 32], 'ClipLimit', 0.01, 'NBins', 64);
       if i<4  
           
           h=figs(:, :, 4);h=h.*mask;h=hist(h(:), 1000);
           dataBrain=histeq(dataBrain, h);         
       else
           %dataBrain=imadjust(dataBrain, [0.05 1-0.05], [0 1]);
       end
       data=dataBrain.*mask+dataSkull.*(~mask);
     %  data=imadjust(data, [0 deltas{i}], [0 1]);
     %  data=imadjust(data, [alphas{i} 1-alphas{i}], [0 1]);
       if i<4
       %    data=resizem(data, [256 256], 'bicubic');
       %    data=resizem(data, [512 512], 'bicubic');
       end
       figs(:, :, i)=data;        
    end    
    h=figure;
    set(h, 'position', [0 100 1440 720]);   
    limits=[80 450 130 330];
    for i=1:4
        subaxis(2, 4, i, 'Margin', 0, 'Spacing', 0.02, 'Padding', 0);       
        colormap(gray);imshow(figs(:, :, i)/1.5);     
        axis off;
        xlim([limits(1) limits(2)]);ylim([limits(3) limits(4)]);
    end
    
    for i=1:3
        subaxis(2, 4, i+4, 'Margin', 0, 'Spacing', 0.02, 'Padding', 0);    
        dif=abs(figs(:, :, i)-figs(:, :, 4));        
        if i==1
            dif=dif*1.7;
        elseif i==2
            dif=dif*1.05*1.7;
        else
            dif=dif*0.95*1.7;
        end
        colormap(gray);imshow(dif);      title(sprintf('norm %d', norm(dif.*mask)));   
        axis off;
        xlim([limits(1) limits(2)]);ylim([limits(3) limits(4)]);
    end

end


function kData=test2(filnam,xdim,ydim)
    info = imfinfo(filnam,'ras');
    headSize=info.FileSize-xdim*ydim*8;
    id=fopen(filnam,'r');
    [header, count] = fread(id, headSize, 'uchar');   
    [kData_1D, count] = fread(id, xdim*ydim*2, 'float');
    fclose(id);
    kReal=kData_1D(1:2:end);
    kImag=kData_1D(2:2:end);
    kData=complex(kReal, kImag);
    kData=reshape(kData,xdim,ydim);
end
%}