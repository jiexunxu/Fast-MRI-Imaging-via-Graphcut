% Compute an initial estimate by finding the analytical solution to
% ||y-Ex||_2^2+lambda*TV(x)_2^2 s.t |Wx|_1<tau
% lambda parameter can be found in the line 
% mu = 0.001;      
function x0=estimateInitialSolution2(data0, sensit0, G, sigma, lambda, maxIter)
    sizeX=[size(data0, 1) size(data0, 2) 4];   
    [nbrList, ~, nbrTermWeights]=computeNbrList(true(sizeX), false(sizeX), [1 1]);
    N=size(nbrList, 1);
    fullNbrMat=sparse([1:N 1:N]', nbrList(:), [nbrTermWeights; -nbrTermWeights], N, sizeX(1)*sizeX(2)*sizeX(3)); 
    fullNbrMatTransp=fullNbrMat';
    conjSensit0=conj(sensit0);    
    wname='db4';decompLvl=3;
    fullNbrMat=lambda*fullNbrMat;      
    data0=reshape(data0, size(data0, 1), size(data0, 2), 4, round(size(data0, 3)/4), size(data0, 4));
    sensit0=reshape(sensit0, size(data0));
    conjSensit0=reshape(conjSensit0, size(data0));
    x0=zeros(size(data0, 1), size(data0, 2), size(data0, 3), size(data0, 4));
    
    parfor iter=1:size(data0, 4)   
        data=squeeze(data0(:, :, :, iter, :));
        sensit=squeeze(sensit0(:, :, :, iter, :));
        conjSensit=squeeze(conjSensit0(:, :, :, iter, :));
        M=numel(data);

        b=zeros(M+size(fullNbrMat, 1), 1);b(1:M)=data(:);
        b=double(b);
            
        W=my_wavedec3(zeros(size(data, 1), size(data, 2), size(data, 3)), decompLvl, wname);
        ylength=0;
        for j=1:length(W.dec)
            ylength=ylength+numel(W.dec{j});
        end
        opts=spgSetParms('iscomplex', 1, 'iterations', maxIter);
        opts.W=W;opts.G=G;opts.sensit=sensit;opts.b=b;opts.fullNbrMat=fullNbrMat;
        opts.conjSensit=conjSensit;opts.fullNbrMatTransp=fullNbrMatTransp;
        opts.decompLvl=decompLvl;opts.wname=wname;opts.M=M;opts.data=data;
        opts.ylength=ylength;
     %   x0Wavelet=spgl1(@HxTV, b, 0, norm(b)*sigma, [], opts);  
        x0Wavelet=my_spgl1(opts, b, 0, norm(b)*sigma, [], opts);  
        xtrunkCounter=1;
        for j=1:length(W.dec)
            sizesRow=max(ceil((j-1)/7), 1);
            sizes=W.sizes(sizesRow, :);
            xtrunk=x0Wavelet(xtrunkCounter:xtrunkCounter+sizes(1)*sizes(2)*sizes(3)-1);
            xtrunkCounter=xtrunkCounter+sizes(1)*sizes(2)*sizes(3);
            W.dec{j}=reshape(xtrunk, sizes);
        end
        x0(:, :, :, iter)=my_waverec3(W);
    end
    x0=reshape(x0, size(x0, 1), size(x0, 2), size(x0, 3)*size(x0, 4));
    
    function y=HxTV(x, mode)
        if mode==1
            % Inverse wavelet transform
            xtrunkCounter=1;
            for i=1:length(W.dec)
                sizesRow=max(ceil((i-1)/7), 1);
                sizes=W.sizes(sizesRow, :);
                xtrunk=x(xtrunkCounter:xtrunkCounter+sizes(1)*sizes(2)*sizes(3)-1);
                xtrunkCounter=xtrunkCounter+sizes(1)*sizes(2)*sizes(3);
                W.dec{i}=reshape(xtrunk, sizes);
            end
            xtmp=my_waverec3(W);
            y=zeros(length(b), 1);
            y(1:M)=Hx(xtmp, G, sensit);
            y(M+1:length(y))=fullNbrMat*xtmp(:);            
        else
            xtmp=Htx(reshape(x(1:M), size(data)), conjSensit, G);
            xtmp=xtmp(:)+fullNbrMatTransp*x(M+1:length(x));
            xtmp=reshape(xtmp, size(data, 1), size(data, 2), size(data, 3));
            W=my_wavedec3(xtmp, decompLvl, wname);
            y=zeros(ylength, 1);
            trunkCounter=1;
            for i=1:length(W.dec)
                y(trunkCounter:trunkCounter+numel(W.dec{i})-1)=W.dec{i}(:);
                trunkCounter=trunkCounter+numel(W.dec{i});
            end
        end
    end

    %{
    % Computes b=H'*data
    b=data;b=Htx(b, sensit, G);b=b(:);
    lambda=lambda*lambda;
    y=sensit;
    x0Size=[size(sensit, 1) size(sensit, 2) size(sensit, 3)];
    numPixels=size(sensit, 1)*size(sensit, 2);
    conjSensit=conj(sensit);
    x0=pcg(@(x) HtHregx(x), b, 1e-6, 100);
    x0=reshape(x0, x0Size);
    
    % Computes y=(H'H+lambda^2*I)x for multiple coils and multiple slices. Both x and y
    % are in image space
    % Note: y is in the IMAGE SPACE, not k-space
    function z=HtHregx(x)
        y=sensit;
        y=bsxfun(@times, y, reshape(x, x0Size));
        y=fft(fft(y, [], 1), [], 2);
        y=bsxfun(@times, y, G);
        y=conj(fft(fft(conj(y), [], 1), [], 2))/numPixels;
        y=y.*conjSensit;
        z=sum(y, 4);
        z=z(:)+lambda*x;
    end
%}
end

function wdec = my_wavedec3(X,level, wname)
%WAVEDEC3 Multilevel 3-D wavelet decomposition.
%   WDEC = WAVEDEC3(X,N,'wname','mode','ExtM') returns the wavelet
%   decomposition of the 3-D array X at level N, using the wavelet 
%   named in string 'wname' (see WFILTERS) or particular wavelet filters
%   you specify, and using a specified DWT extension mode (see DWTMODE).
%   WDEC = WAVEDEC3(X,N,'wname') uses the default extension mode: 'sym'.
%   
%   N must be a strictly positive integer (see WMAXLEV).
%
%   WDEC is the output decomposition structure, with the following fields:
%     sizeINI: contains the size of the 3-D array X.
%     level:   contains the level of the decomposition.
%     mode:    contains the name of the wavelet transform extension mode.
%     filters: is a structure with 4 fields LoD, HiD, LoR, HiR which
%              contain the filters used for DWT.
%     dec:     is a Nx1 cell array containing the coefficients 
%              of the decomposition. N is equal to 7*level+1.
%              dec(1) contains the low pass component (approximation),
%              dec(k+2),...,dec(k+8) contain the components of level L,
%              L = (level-k) with k = 0,...,level-1, in the following
%              order: 'LLH','LHL','LHH','HLL','HLH','HHL','HHH'.
%     sizes:   contains the successive sizes of the decomposition
%              components.
%
%   Examples:
%       M = magic(8);
%       X = repmat(M,[1 1 8]);
%       wd1 = wavedec3(X,1,'db1')
%       [LoD,HiD,LoR,HiR] = wfilters('db2');
%       wd2 = wavedec3(X,2,{LoD,HiD,LoR,HiR})
%       wd3 = wavedec3(X,2,{LoD,HiD,LoR,HiR},'mode','per')
%
%   See also dwtmode, dwt3, waverec3, waveinfo.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 09-Dec-2008.
%   Last Revision: 04-Oct-2009.
%   Copyright 1995-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $ $Date: 2009/11/13 05:31:34 $

% Check arguments.
LD=[-1.05974017849973e-002 3.28830116669829e-002 3.08413818359870e-002 -1.87034811718881e-001 -2.79837694169838e-002 6.30880767929590e-001 7.14846570552542e-001 2.30377813308855e-001];
HD=[-2.30377813308855e-001 7.14846570552542e-001 -6.30880767929590e-001 -2.79837694169838e-002 1.87034811718881e-001 3.08413818359870e-002 -3.28830116669829e-002 -1.05974017849973e-002];
LR=[2.30377813308855e-001 7.14846570552542e-001 6.30880767929590e-001 -2.79837694169838e-002 -1.87034811718881e-001 3.08413818359870e-002 3.28830116669829e-002 -1.05974017849973e-002];
HR=[-1.05974017849973e-002 -3.28830116669829e-002 3.08413818359870e-002 1.87034811718881e-001 -2.79837694169838e-002 -6.30880767929590e-001 7.14846570552542e-001 -2.30377813308855e-001];

% Initialization.
if isempty(X) , wdec = {}; return; end
sizes = zeros(level+1,3);
sizes(level+1,1:3) = size(X);
for k=1:level
    wdec = my_dwt3(X, LD,HD,LR,HR);
    X = wdec.dec{1,1,1};
    if length(size(X))>2
        sizes(level+1-k,1:3) = size(X);
    else
        sizes(level+1-k,1:3) = ceil(sizes(level+2-k,1:3)/2);
    end
    wdec.dec = reshape(wdec.dec,8,1,1);
    if k>1
        cfs(1) = [];
        cfs = cat(1,wdec.dec,cfs);
    else
        cfs = wdec.dec;
    end
end
wdec.sizeINI = sizes(end,:);
wdec.level = level;
wdec.dec   = cfs;
wdec.sizes = sizes;
wdec = orderfields(wdec,{'sizeINI','level','filters','mode','dec','sizes'});
end
function wt = my_dwt3(X, LD, HD, LR, HR)
%DWT3 Single-level discrete 3-D wavelet transform.
%   DWT3 performs a single-level 3-D wavelet decomposition
%   with respect to either a particular wavelet ('wname',
%   see WFILTERS for more information) or particular wavelet 
%   decomposition and reconstruction filters you specify, and 
%   using a specified DWT extension mode (see DWTMODE).
%
%   WT = DWT3(X,'wname','mode','ExtM') returns the 3-D wavelet transform
%   of the 3-D array X, 'wname' is a string containing the wavelet 
%   name and 'ExtM' is a string containing the extension mode.
%   WT = DWT3(X,'wname') uses the default extension mode: 'sym'.
%   
%   WT is a structure with the following fields:
%     sizeINI: contains the size of the 3-D array X.
%     mode:    contains the name of the wavelet transform extension mode.
%     filters: is a structure with 4 fields LoD, HiD, LoR, HiR which
%              contain the filters used for DWT.
%     dec:     is a 2x2x2 cell array containing the coefficients 
%              of the decomposition.
%              dec{i,j,k} for i,j,k = 1 or 2 contains the coefficients
%              obtained by low-pass filtering (for i or j or k = 1)
%              or high-pass filtering (for i or j or k = 2).              
%
%   Instead of a single wavelet, you may specify three wavelets (i.e. 
%   one wavelet for each direction):
%   WT = DWT3(X,W,...) with W = {'wname1','wname2','wname3'} or W a 
%   structure with 3 fields 'w1', 'w2', 'w3' containing strings which
%   are the names of wavelets.
%
%   Instead of wavelets you may specify filters: 4 filters (2 for
%   decomposition and 2 for reconstruction) or 3x4 filters (one 
%   quadruplet by direction): WT = DWT3(X,WF,...)
%   Where WF must be a cell array (1x4) or (3x4) : {LoD,HiD,LoR,HiR},
%   or a structure with the four fields 'LoD', 'HiD', 'LoR', 'HiR'.
%
%   Examples:
%       X = reshape(1:64,4,4,4)
%       wt = dwt3(X,'db1')
%       [LoD,HiD,LoR,HiR] = wfilters('db2');
%       wt = dwt3(X,{LoD,HiD,LoR,HiR})
%       WS = struct('w1','db1','w2','db2','w3','db1');
%       wt = dwt3(X,WS,'mode','per')
%       WF = wt.filters;
%       wtBIS = dwt3(X,WF,'mode','sym')
%
%   See also dwtmode, idwt3, wavedec3, waverec3, waveinfo.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 09-Dec-2008.
%   Last Revision: 21-Oct-2009.
%   Copyright 1995-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $ $Date: 2009/11/13 05:31:07 $

% Check arguments.

LoD = cell(1,3); HiD = cell(1,3); LoR = cell(1,3); HiR = cell(1,3);    
for k = 1:3
    LoD{k} = LD; HiD{k} = HD; LoR{k} = LR; HiR{k} = HR;
end
    
sX = size(X);

% Check arguments for Extension.
dwtEXTM = 'sym';

X = double(X);
dec = cell(2,2,2);
permVect = [];
[a_Lo,d_Hi] = wdec1D(X,LoD{1},HiD{1},permVect,dwtEXTM);
permVect = [2,1,3];
[aa_Lo_Lo,da_Lo_Hi] = wdec1D(a_Lo,LoD{2},HiD{2},permVect,dwtEXTM);
[ad_Hi_Lo,dd_Hi_Hi] = wdec1D(d_Hi,LoD{2},HiD{2},permVect,dwtEXTM);
permVect = [1,3,2];
[dec{1,1,1},dec{1,1,2}] = wdec1D(aa_Lo_Lo,LoD{3},HiD{3},permVect,dwtEXTM);
[dec{2,1,1},dec{2,1,2}] = wdec1D(da_Lo_Hi,LoD{3},HiD{3},permVect,dwtEXTM);
[dec{1,2,1},dec{1,2,2}] = wdec1D(ad_Hi_Lo,LoD{3},HiD{3},permVect,dwtEXTM);
[dec{2,2,1},dec{2,2,2}] = wdec1D(dd_Hi_Hi,LoD{3},HiD{3},permVect,dwtEXTM);
wt.sizeINI = sX;
wt.filters.LoD = LoD;
wt.filters.HiD = HiD;
wt.filters.LoR = LoR;
wt.filters.HiR = HiR;
wt.mode = dwtEXTM;
wt.dec = dec;
end
%-----------------------------------------------------------------------%
function [L,H] = wdec1D(X,Lo,Hi,perm,dwtEXTM)

if ~isempty(perm) , X = permute(X,perm); end
sX = size(X);
if length(sX)<3 , sX(3) = 1; end

lf = length(Lo);
lx = sX(2);
lc = lx+lf-1;
if lx<lf+1
    nbAdd = lf-lx+1;
    switch dwtEXTM
        case {'sym','symh','symw','asym','asymh','asymw','ppd'}
            Add = zeros(sX(1),nbAdd,sX(3));
            X = [Add , X , Add];
    end
end

switch dwtEXTM
    case 'zpd'             % Zero extension.
        
    case {'sym','symh'}    % Symmetric extension (half-point).
        X = [X(:,lf-1:-1:1,:) , X , X(:,end:-1:end-lf+1,:)];
        
    case 'sp0'             % Smooth extension of order 0.
        X = [X(:,ones(1,lf-1),:) , X , X(:,lx*ones(1,lf-1),:)];
        
    case {'sp1','spd'}     % Smooth extension of order 1.
        Z = zeros(sX(1),sX(2)+ 2*lf-2,sX(3));
        Z(:,lf:lf+lx-1,:) = X;
        last = sX(2)+lf-1;
        for k = 1:lf-1
            Z(:,last+k,:) = 2*Z(:,last+k-1,:)- Z(:,last+k-2,:);
            Z(:,lf-k,:)   = 2*Z(:,lf-k+1,:)- Z(:,lf-k+2,:);
        end
        X = Z; clear Z;
        
    case 'symw'            % Symmetric extension (whole-point).
        X = [X(:,lf:-1:2,:) , X , X(:,end-1:-1:end-lf,:)];
        
    case {'asym','asymh'}  % Antisymmetric extension (half-point).
        X = [-X(:,lf-1:-1:1,:) , X , -X(:,end:-1:end-lf+1,:)];        
        
    case 'asymw'           % Antisymmetric extension (whole-point).
        X = [-X(:,lf:-1:2,:) , X , -X(:,end-1:-1:end-lf,:)];

    case 'rndu'            % Uniformly randomized extension.
        X = [randn(sX(1),lf-1,sX(3)) , X , randn(sX(1),lf-1,sX(3))];        
                        
    case 'rndn'            % Normally randomized extension.
        X = [randn(sX(1),lf-1,sX(3)) , X , randn(sX(1),lf-1,sX(3))];        
                
    case 'ppd'             % Periodized extension (1).
        X = [X(:,end-lf+2:end,:) , X , X(:,1:lf-1,:)];
        
    case 'per'             % Periodized extension (2).
        if rem(lx,2) , X = [X , X(:,end,:)]; lx = lx + 1; end
        I = [lx-lf+1:lx , 1:lx , 1:lf];
        if lx<lf
            I = mod(I,lx);
            I(I==0) = lx;
        end
        X = X(:,I,:);
end
L = convn(X,Lo);
H = convn(X,Hi);
clear X
switch dwtEXTM
    case 'zpd'
    otherwise
        lenL = size(L,2);
        first = lf; last = lenL-lf+1;
        L = L(:,first:last,:); H = H(:,first:last,:);
        lenL = size(L,2);
        first = 1+floor((lenL-lc)/2);  last = first+lc-1;
        L = L(:,first:last,:); H = H(:,first:last,:);
end
L = L(:,2:2:end,:);
H = H(:,2:2:end,:);
if isequal(dwtEXTM,'per')
    last = ceil(lx/2);
    L = L(:,1:last,:);
    H = H(:,1:last,:);
end

if ~isempty(perm)
    L = permute(L,perm);
    H = permute(H,perm);
end
%-----------------------------------------------------------------------%
end

function X = my_waverec3(wdec)
%WAVEREC3 Multilevel 3-D wavelet reconstruction.
%   WAVEREC3 performs a multilevel 3-D wavelet reconstruction
%   starting from a multilevel 3-D wavelet decomposition.
%
%   X = WAVEREC3(WDEC) reconstructs the 3-D array X based
%   on the multilevel wavelet decomposition structure WDEC.
%
%   In addition, you can use WAVEREC3 to simply extract coefficients 
%   from a 3-D wavelet decomposition (see below).
%
%   WDEC is a structure with the following fields:
%     sizeINI: contains the size of the 3-D array X.
%     level:   contains the level of the decomposition.
%     mode:    contains the name of the wavelet transform extension mode.
%     filters: is a structure with 4 fields LoD, HiD, LoR, HiR which
%              contain the filters used for DWT.
%     dec:     is a Nx1 cell array containing the coefficients 
%              of the decomposition. N is equal to 7*level+1.
%              dec(1) contains the low pass component (approximation),
%              dec(k+2),...,dec(k+8) contain the components of level L,
%              L = (level-k) with k = 0,...,level-1, in the following
%              order: 'LLH','LHL','LHH','HLL','HLH','HHL','HHH'.
%     sizes:   contains the successive sizes of the decomposition
%              components.
%
%   C = WAVEREC3(WDEC,TYPE,N) reconstructs the multilevel
%   components at level N of a 3-D wavelet decomposition. N must be
%   a positive integer less or equal to the level of the decomposition.
%   The valid values for TYPE are:
%       - A group of 3 chars 'xyz', one per direction, with 'x','y' and 'z' 
%         in the set {'a','d','l','h'} or in the corresponding upper case  
%         set {'A','D','L','H'}), where 'A' (or 'L') stands for low pass 
%         filter and 'D' (or 'H') stands for high pass filter.
%       - The char 'd' (or 'h' or 'D' or 'H') gives directly the sum of 
%         all the components different from the low pass one.
%       - The char 'a' (or 'l' or 'A' or 'L') gives the low pass 
%         component (the approximation at level N).
%
%   For extraction purpose, the valid values for TYPE are the same as
%   above prefixed by 'c' or 'C'.
%
%   X = WAVEREC3(WDEC,'a',0) or X = WAVEREC3(WDEC,'ca',0) is equivalent
%   to X = WAVEREC3(WDEC).
%
%	C = WAVEREC3(WDEC,TYPE) is equivalent to C = WAVEREC3(WDEC,TYPE,N)
%   with N equal to the level of the decomposition.
%
%   Examples:
%       M = magic(8);
%       X = repmat(M,[1 1 8]);
%       wd = wavedec3(X,2,'db2','mode','per');
%       XR = waverec3(wd);
%       err1 = max(abs(X(:)-XR(:)))
%       A = waverec3(wd,'aaa');
%       CA = waverec3(wd,'ca');
%       D = waverec3(wd,'d');
%       err2 = max(abs(X(:)-A(:)-D(:)))
%       A1 = waverec3(wd,'aaa',1);
%       D1 = waverec3(wd,'d',1);
%       err3 = max(abs(X(:)-A1(:)-D1(:)))
%       DDD = waverec3(wd,'ddd',2);
%       disp(DDD(:,:,1))
%
%   See also idwt3, wavedec3, waveinfo.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 09-Dec-2008.
%   Last Revision: 12-Oct-2009.
%   Copyright 1995-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $ $Date: 2009/11/13 05:31:36 $

% Check arguments.

% Initialization.
cfs   = wdec.dec;
sizes = wdec.sizes;
level = wdec.level;
nbREC_UP = level;

idxBeg = 1;
for k=1:nbREC_UP
    idxEnd = idxBeg+7;
    wdec.dec = reshape(cfs(idxBeg:idxEnd),2,2,2);
    X = my_idwt3(wdec,sizes(k+1,:));
    cfs{idxEnd} = X;
    idxBeg = idxEnd;
end
end

function X = my_idwt3(wt,varargin)
%IDWT3 Single-level inverse discrete 3-D wavelet transform.
%   IDWT3 performs a single-level 3-D wavelet reconstruction
%   starting from a single-level 3-D wavelet decomposition.
%
%   X = IDWT3(WT) computes the single-level reconstructed 3-D array
%   X based on 3-D wavelet decomposition contained in the structure
%   WT which contains the following fields:
%     sizeINI: contains the size of the 3-D array X.
%     mode:    contains the name of the wavelet transform extension mode.
%     filters: is a structure with 4 fields LoD, HiD, LoR, HiR which
%              contain the filters used for DWT.
%     dec:     is a 2x2x2 cell array containing the coefficients 
%              of the decomposition.
%              dec{i,j,k} , i,j,k = 1 or 2 contains the coefficients
%              obtained by low-pass filtering (for i or j or k = 1)
%              or high-pass filtering (for i or j or k = 2).              
%
%   C = IDWT3(WT,TYPE) allows to compute the single-level reconstructed 
%   component based on the 3-D wavelet decomposition. 
%   The valid values for TYPE are:
%       - A group of 3 chars 'xyz', one per direction, with 'x','y' and 'z' 
%         in the set {'a','d','l','h'} or in the corresponding upper case  
%         set {'A','D','L','H'}), where 'A' (or 'L') stands for low pass 
%         filter and 'D' (or 'H') stands for high pass filter.
%       - The char 'd' (or 'h' or 'D' or 'H') gives directly the sum of 
%         all the components different from the low pass one.
%
%   Examples:
%       X  = reshape(1:64,4,4,4);
%       wt = dwt3(X,'db1');
%       XR = idwt3(wt);
%       A  = idwt3(wt,'aaa');
%       D  = idwt3(wt,'d');
%       ADA  = idwt3(wt,'ada');
%
%   See also dwt3, wavedec3, waverec3.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 09-Dec-2008.
%   Last Revision: 21-Oct-2009.
%   Copyright 1995-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $ $Date: 2009/11/13 05:31:13 $

% Check arguments.
nbIn = nargin;
msg = nargchk(1,2,nbIn); %#ok<NCHK>

s   = wt.sizeINI;
dec = wt.dec;
Lo  = wt.filters.LoR;
Hi  = wt.filters.HiR;
dwtEXTM = wt.mode;
perFLAG = isequal(dwtEXTM,'per');

lf = zeros(1,3);
for k = 1:3 , lf(k) = length(Lo{k}); end

k = 1;
while k<=length(varargin)
        s = varargin{k};
        k = k+1;
end

% Reconstruction.
perm = [1,3,2];
V = cell(2,2);
for i = 1:2    
    for j = 1:2
        V{j,i} = wrec1D(dec{i,j,1},Lo{3},perm,perFLAG,s) + ...
                 wrec1D(dec{i,j,2},Hi{3},perm,perFLAG,s);
    end
end
perm = [2,1,3];
W = cell(1,2);
for i = 1:2
    W{i} = wrec1D(V{i,1},Lo{2},perm,perFLAG,s) + ...
        wrec1D(V{i,2},Hi{2},perm,perFLAG,s);
end

% Last reconstruction.
X = wrec1D(W{1},Lo{1},[],perFLAG,s) + wrec1D(W{2},Hi{1},[],perFLAG,s);
end
%-----------------------------------------------------------------------%
function X = wrec1D(X,F,perm,perFLAG,s)

if ~isempty(perm)
    X = permute(X,perm);
    s = s(perm);
end
if perFLAG
    lf = length(F);
    lx = size(X,2);
    nb = fix(lf/2-1);
    idxAdd = 1:nb;
    if nb>lx
        idxAdd = rem(idxAdd,lx);
        idxAdd(idxAdd==0) = lx;
    end
    X = [X X(:,idxAdd,:)];
end
sX = size(X);
if length(sX)<3 , sX(3) = 1; end
Z = zeros(sX(1),2*sX(2)-1,sX(3));
Z(:,1:2:end,:) = X;
X = convn(Z,F);

sX = size(X,2);
F  = floor((sX-s)/2);
C  = ceil((sX-s)/2);
X  = X(:,1+F(2):end-C(2),:);

if ~isempty(perm) , X = permute(X,perm); end
%-----------------------------------------------------------------------%
end
