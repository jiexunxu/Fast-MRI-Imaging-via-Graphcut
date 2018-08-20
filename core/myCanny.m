% reference http://www.ncbi.nlm.nih.gov/pubmed/19965004
%{
function e = myCanny(varargin)
    e = myCannyPrivate(varargin);
    %{
    img=varargin{1};
    e=zeros(size(img, 1), size(img, 2), size(img, 3));
    for i=1:size(img, 3)
        varargin{1}=img(:, :, i);
        e1 = myCannyPrivate(varargin);
        e(:, :, i) = e1;
    end
    %}
%}

function e = myCanny(varargin)
[a,method,thresh,sigma,thinning,H,kx,ky] = parse_inputs(varargin{:});

% Check that the user specified a valid number of output arguments
if ~any(strcmp(method,{'sobel','roberts','prewitt'})) && (nargout>2)
    eid = sprintf('Images:%s:tooManyOutputs', mfilename);
    msg = 'Too many output arguments';
    error(eid,'%s',msg);
end

% Transform to a double precision intensity image if necessary
if ~isa(a,'double') && ~isa(a,'single')
    a = im2single(a);
end

[m,n] = size(a);

% The output edge map:
e = false(m,n);

if strcmp(method,'canny')
    % Magic numbers
    PercentOfPixelsNotEdges = .7; % Used for selecting thresholds
    ThresholdRatio = .3;          % Low thresh is this fraction of the high.
    
    % Calculate gradients using a derivative of Gaussian filter
    [dx, dy] = smoothGradient(a, sigma);
    
    % Calculate Magnitude of Gradient
    
    theta=atan(dy./dx);
    magGrad=zeros(size(dy, 1), size(dy, 2), 2);
    magGrad(:, :, 1)=abs(cos(theta).*dx);
    magGrad(:, :, 2)=abs(sin(theta).*dy);
    magGrad=max(magGrad, [], 3);  
    magGrad2 = hypot(dx, dy);
    
    % Normalize for threshold selection
    magmax = max(magGrad(:));magmax2=max(magGrad2(:));
    magGrad = magGrad / magmax;
    magGrad2=magGrad2/magmax2;
    
    % Determine Hysteresis Thresholds
    [lowThresh, highThresh] = selectThresholds(thresh, magGrad, PercentOfPixelsNotEdges, ThresholdRatio, mfilename);
    e1 = thinAndThreshold(e, dx, dy, magGrad, lowThresh, highThresh);
    
    [lowThresh, highThresh] = selectThresholds(thresh, magGrad2, PercentOfPixelsNotEdges, ThresholdRatio, mfilename);
    e2 = thinAndThreshold(e, dx, dy, magGrad2, lowThresh, highThresh);
    
    e=e1+e2; e(e==2)=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : cannyFindLocalMaxima
%
function idxLocalMax = cannyFindLocalMaxima(direction,ix,iy,mag)
%
% This sub-function helps with the non-maximum suppression in the Canny
% edge detector.  The input parameters are:
%
%   direction - the index of which direction the gradient is pointing,
%               read from the diagram below. direction is 1, 2, 3, or 4.
%   ix        - input image filtered by derivative of gaussian along x
%   iy        - input image filtered by derivative of gaussian along y
%   mag       - the gradient magnitude image
%
%    there are 4 cases:
%
%                         The X marks the pixel in question, and each
%         3     2         of the quadrants for the gradient vector
%       O----0----0       fall into two cases, divided by the 45
%     4 |         | 1     degree line.  In one case the gradient
%       |         |       vector is more horizontal, and in the other
%       O    X    O       it is more vertical.  There are eight
%       |         |       divisions, but for the non-maximum suppression
%    (1)|         |(4)    we are only worried about 4 of them since we
%       O----O----O       use symmetric points about the center pixel.
%        (2)   (3)


[m,n] = size(mag);

% Find the indices of all points whose gradient (specified by the
% vector (ix,iy)) is going in the direction we're looking at.

switch direction
    case 1
        idx = find((iy<=0 & ix>-iy)  | (iy>=0 & ix<-iy));
    case 2
        idx = find((ix>0 & -iy>=ix)  | (ix<0 & -iy<=ix));
    case 3
        idx = find((ix<=0 & ix>iy) | (ix>=0 & ix<iy));
    case 4
        idx = find((iy<0 & ix<=iy) | (iy>0 & ix>=iy));
end

% Exclude the exterior pixels
if ~isempty(idx)
    v = mod(idx,m);
    extIdx = (v==1 | v==0 | idx<=m | (idx>(n-1)*m));
    idx(extIdx) = [];
end

ixv = ix(idx);
iyv = iy(idx);
gradmag = mag(idx);

% Do the linear interpolations for the interior pixels
switch direction
    case 1
        d = abs(iyv./ixv);
        gradmag1 = mag(idx+m).*(1-d) + mag(idx+m-1).*d;
        gradmag2 = mag(idx-m).*(1-d) + mag(idx-m+1).*d;
    case 2
        d = abs(ixv./iyv);
        gradmag1 = mag(idx-1).*(1-d) + mag(idx+m-1).*d;
        gradmag2 = mag(idx+1).*(1-d) + mag(idx-m+1).*d;
    case 3
        d = abs(ixv./iyv);
        gradmag1 = mag(idx-1).*(1-d) + mag(idx-m-1).*d;
        gradmag2 = mag(idx+1).*(1-d) + mag(idx+m+1).*d;
    case 4
        d = abs(iyv./ixv);
        gradmag1 = mag(idx-m).*(1-d) + mag(idx-m-1).*d;
        gradmag2 = mag(idx+m).*(1-d) + mag(idx+m+1).*d;
end
idxLocalMax = idx(gradmag>=gradmag1 & gradmag>=gradmag2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : parse_inputs
%
function [I,Method,Thresh,Sigma,Thinning,H,kx,ky] = parse_inputs(varargin)
% OUTPUTS:
%   I      Image Data
%   Method Edge detection method
%   Thresh Threshold value
%   Sigma  standard deviation of Gaussian
%   H      Filter for Zero-crossing detection
%   kx,ky  From Directionality vector

error(nargchk(1,5,nargin,'struct'));

I = varargin{1};

% iptcheckinput(I,{'numeric','logical'},{'nonsparse','2d'},mfilename,'I',1);

% Defaults
Method='sobel';
Thresh=[];
Direction='both';
Thinning=true;
Sigma=2;
H=[];
K=[1 1];

methods = {'canny','canny_old','prewitt','sobel','marr-hildreth','log','roberts','zerocross'};
directions = {'both','horizontal','vertical'};
options = {'thinning','nothinning'};

% Now parse the nargin-1 remaining input arguments

% First get the strings - we do this because the interpretation of the
% rest of the arguments will depend on the method.
nonstr = [];   % ordered indices of non-string arguments
for i = 2:nargin
    if ischar(varargin{i})
        str = lower(varargin{i});
        j = find(strcmp(str,methods));
        k = find(strcmp(str,directions));
        l = find(strcmp(str,options));
        if ~isempty(j)
            Method = methods{j(1)};
            if strcmp(Method,'marr-hildreth')
                wid = sprintf('Images:%s:obsoleteMarrHildrethSyntax', mfilename);
                msg = '''Marr-Hildreth'' is an obsolete syntax, use ''LoG'' instead.';
                warning(wid,'%s',msg);
            end
        elseif ~isempty(k)
            Direction = directions{k(1)};
        elseif ~isempty(l)
            if strcmp(options{l(1)},'thinning')
                Thinning = true;
            else
                Thinning = false;
            end
        else
            eid = sprintf('Images:%s:invalidInputString', mfilename);
            msg = sprintf('%s%s%s', 'Invalid input string: ''', varargin{i},'''.');
            error(eid,'%s',msg);
        end
    else
        nonstr = [nonstr i]; %#ok<AGROW>
    end
end

% Now get the rest of the arguments

eid_invalidArgs = sprintf('Images:%s:invalidInputArguments', mfilename);
msg_invalidArgs = 'Invalid input arguments';

switch Method
    
    case {'prewitt','sobel','roberts'}
        threshSpecified = 0;  % Threshold is not yet specified
        for i = nonstr
            if numel(varargin{i})<=1 && ~threshSpecified % Scalar or empty
                Thresh = varargin{i};
                threshSpecified = 1;
            elseif numel(varargin{i})==2  % The dreaded K vector
                wid = sprintf('Images:%s:obsoleteKDirectionSyntax', mfilename);
                msg = sprintf('%s%s%s', 'BW = EDGE(... , K) is an obsolete syntax. ',...
                    'Use BW = EDGE(... , DIRECTION),',...
                    ' where DIRECTION is a string.');
                warning(wid,'%s',msg);
                K=varargin{i};
            else
                error(eid_invalidArgs,msg_invalidArgs);
            end
        end
        
    case {'canny','canny_old'}
        % Default Std dev of gaussian for canny
        if (strcmp(Method, 'canny'))
            % The canny_old method smooths the input image twice
            % Use sqrt(2) to achieve similar results to the canny_old
            % method
            Sigma = sqrt(2);  
        else
            Sigma = 1.0;
        end
        threshSpecified = 0;  % Threshold is not yet specified
        for i = nonstr
            if numel(varargin{i})==2 && ~threshSpecified
                Thresh = varargin{i};
                threshSpecified = 1;
            elseif numel(varargin{i})==1
                if ~threshSpecified
                    Thresh = varargin{i};
                    threshSpecified = 1;
                else
                    Sigma = varargin{i};
                end
            elseif isempty(varargin{i}) && ~threshSpecified
                threshSpecified = 1;
            else
                error(eid_invalidArgs,msg_invalidArgs);
            end
        end
        
    case 'log'
        threshSpecified = 0;  % Threshold is not yet specified
        for i = nonstr
            if numel(varargin{i})<=1  % Scalar or empty
                if ~threshSpecified
                    Thresh = varargin{i};
                    threshSpecified = 1;
                else
                    Sigma = varargin{i};
                end
            else
                error(eid_invalidArgs,msg_invalidArgs);
            end
        end
        
    case 'zerocross'
        threshSpecified = 0;  % Threshold is not yet specified
        for i = nonstr
            if numel(varargin{i})<=1 && ~threshSpecified % Scalar or empty
                Thresh = varargin{i};
                threshSpecified = 1;
            elseif numel(varargin{i}) > 1 % The filter for zerocross
                H = varargin{i};
            else
                error(eid_invalidArgs,msg_invalidArgs);
            end
        end
        
    case 'marr-hildreth'
        for i = nonstr
            if numel(varargin{i})<=1  % Scalar or empty
                Thresh = varargin{i};
            elseif numel(varargin{i})==2  % The dreaded K vector
                wid = sprintf('Images:%s:dirFactorHasNoEffectOnMarrHildreth', mfilename);
                msg = 'The [kx ky] direction factor has no effect for ''Marr-Hildreth''.';
                warning(wid,'%s',msg);
            elseif numel(varargin{i}) > 2 % The filter for zerocross
                H = varargin{i};
            else
                error(eid_invalidArgs,msg_invalidArgs);
            end
        end
        
    otherwise
        error(eid_invalidArgs,msg_invalidArgs);
        
end

if Sigma<=0
    eid = sprintf('Images:%s:sigmaMustBePositive', mfilename);
    msg = 'Sigma must be positive';
    error(eid,'%s',msg);
end

switch Direction
    case 'both',
        kx = K(1); ky = K(2);
    case 'horizontal',
        kx = 0; ky = 1; % Directionality factor
    case 'vertical',
        kx = 1; ky = 0; % Directionality factor
    otherwise
        eid = sprintf('Images:%s:badDirectionString', mfilename);
        msg = 'Unrecognized direction string';
        error(eid,'%s',msg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : smoothGradient
%
function [GX, GY] = smoothGradient(I, sigma)

% Create an even-length 1-D separable Derivative of Gaussian filter
xKernel=[1 sqrt(2) 2 0 -2 -sqrt(2) -1; ...
         sqrt(2) 2 2*sqrt(2) 0 -2*sqrt(2) -2 -sqrt(2); ...
         2 2*sqrt(2) 4 0 -4 -2*sqrt(2) -2;...
         sqrt(2) 2 2*sqrt(2) 0 -2*sqrt(2) -2 -sqrt(2); ...
         1 sqrt(2) 2 0 -2 -sqrt(2) -1; ];
yKernel=xKernel';

% Compute smoothed numerical gradient of image I along x (horizontal)
% direction. GX corresponds to dG/dx, where G is the Gaussian Smoothed
% version of image I.
GX = imfilter(I, xKernel, 'conv', 'replicate');

% Compute smoothed numerical gradient of image I along y (vertical)
% direction. GY corresponds to dG/dy, where G is the Gaussian Smoothed
% version of image I.
GY = imfilter(I, yKernel, 'conv', 'replicate');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : selectThresholds
%
function [lowThresh, highThresh] = selectThresholds(thresh, magGrad, PercentOfPixelsNotEdges, ThresholdRatio, mfilename)

[m,n] = size(magGrad);

% Select the thresholds
if isempty(thresh)
    counts=imhist(magGrad, 64);
    highThresh = find(cumsum(counts) > PercentOfPixelsNotEdges*m*n,...
        1,'first') / 64;
    lowThresh = ThresholdRatio*highThresh;
elseif length(thresh)==1
    highThresh = thresh;
    if thresh>=1
        eid = sprintf('Images:%s:thresholdMustBeLessThanOne', mfilename);
        msg = 'The threshold must be less than 1.';
        error(eid,'%s',msg);
    end
    lowThresh = ThresholdRatio*thresh;
elseif length(thresh)==2
    lowThresh = thresh(1);
    highThresh = thresh(2);
    if (lowThresh >= highThresh) || (highThresh >= 1)
        eid = sprintf('Images:%s:thresholdOutOfRange', mfilename);
        msg = 'Thresh must be [low high], where low < high < 1.';
        error(eid,'%s',msg);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : thinAndThreshold
%
function H = thinAndThreshold(E, dx, dy, magGrad, lowThresh, highThresh)

% Perform Non-Maximum Suppression Thining and Hysteresis Thresholding of Edge
% Strength
    
% We will accrue indices which specify ON pixels in strong edgemap
% The array e will become the weak edge map.
idxStrong = [];
for dir = 1:4
    idxLocalMax = cannyFindLocalMaxima(dir,dx,dy,magGrad);
    idxWeak = idxLocalMax(magGrad(idxLocalMax) > lowThresh);
    E(idxWeak)=1;
    idxStrong = [idxStrong; idxWeak(magGrad(idxWeak) > highThresh)]; %#ok<AGROW>
end

[m,n] = size(E);

if ~isempty(idxStrong) % result is all zeros if idxStrong is empty
    rstrong = rem(idxStrong-1, m)+1;
    cstrong = floor((idxStrong-1)/m)+1;
    H = bwselect(E, cstrong, rstrong, 8);
else
    H = zeros(m, n);
end