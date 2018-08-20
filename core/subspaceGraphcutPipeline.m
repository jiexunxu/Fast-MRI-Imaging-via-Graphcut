function [x, x0]=subspaceGraphcutPipeline(data, sensitivities, G, varargin)
    x0=[];voxelMask=[];
    if nargin==5
        x0=varargin{1};voxelMask=varargin{2};
    elseif nargin==4
        x0=varargin{1};
    end
    params=subspaceGraphcutParameters();
    
    [x0, voxelMask, segVoxelIndices]=graphcutInit(data, sensitivities, G, x0, voxelMask, params);
    x=subspaceGraphcutPipelineHelper(data, sensitivities, G, x0, voxelMask, segVoxelIndices, params);    
end

function [x0, voxelMask, segVoxelIndices]=graphcutInit(data, sensitivities, G, x0, voxelMask, params)
    if isempty(x0)   
        if strcmpi(params.initializationMethod, 'SENSE')
            x0=computeSENSE(data, sensitivities, G, params.SENSEweight, '');
        elseif strcmpi(params.initializationMethod, 'SENSE_TV')
            x0=computeSENSE(double(data), double(sensitivities), G, params.SENSEweight, 'tv');
        elseif strcmpi(params.initializationMethod, 'CS')
            x0=estimateInitialSolution(data, sensitivities, G, params.L1Sigma, params.L1SmoothnessWeight, params.L1MaxIter);       
        end
    end
    if isempty(voxelMask)
        if strcmpi(params.segmentationMethod, 'SPM8Volumetrics')
            [~, ~, segMap]=spmVolumetrics(x0);
            voxelMask=(segMap>0);
            segVoxelIndices=refineSPMSegMap(segMap);
        else
           % voxelMask=generateVoxelMask(x0);
           voxelMask=true(size(x0));
           segVoxelIndices=segmentImageAndMat(abs(x0), voxelMask, params.graphSegmentationMinSizeRatio);             
        end
    else
        segVoxelIndices=segmentImageAndMat(abs(x0), voxelMask, params.graphSegmentationMinSizeRatio); 
    end
    x0=abs(reshape(double(x0), size(data, 1), size(data, 2), size(data, 3)));      
    fftw('planner', 'exhaustive');    
end

function x=subspaceGraphcutPipelineHelper(data, sensitivities, G, x0, voxelMask, segVoxelIndices, params)    
    objParams=[params.dataTermWeight params.EdEsRatio params.truncationFactor params.smoothnessPower];
    
    edgeMap=constructEdgeMap(x0, params.edgeMapThreshold);
    edgeMap=false(size(x0));
    [nbrList, ~, nbrTermWeights]=computeNbrList(voxelMask, edgeMap, [min(objParams(3), 1) 1]);
    N=size(nbrList, 1);
    fullNbrMat=sparse([1:N 1:N]', nbrList(:), [nbrTermWeights; -nbrTermWeights], N, numel(voxelMask));
%{
    xUnsharp=unsharpEnhancement(x0, voxelMask); 
    disp('Finish unsharp masking');
    V=generateWaveletFilterBank(xUnsharp, decomp_level, wname, voxelMask, segVoxelIndices);
    disp('Finish wavelet subspace');
    
    xWavelet=waveletDenoising(V, xUnsharp, sensitivities, G, data, objParams, nbrTermWeights, fullNbrMat); 

    disp('Finish wavelet denoising');
    clear V;
%}
 %   xWavelet=unsharpEnhancement(x0, voxelMask);
    V=generateHighFreqSubspaceVector(x0, voxelMask, segVoxelIndices, params.BFsimgas, ...
        params.BFsigmaWeight, params.use3DBF, G); 
    disp('Finish bilateral subspace');
    VEEV=computeUnaryAndBinaryDataTerms(V, sensitivities, G);
    disp('Finish VEEV');
    opts=struct('inputMatData', struct('data', data, 'sensit', sensitivities, 'G', G), ...
        'subspaceMove', struct('x0', x0, 'subspace', V, 'iterCount', params.graphcutIterCount, ...
        'objParams', objParams, 'dataTermMatrix', VEEV, 'nbrList', nbrList, 'nbrTermWeights', nbrTermWeights, ...
        'fullNbrMat', fullNbrMat, 'numSegs', length(segVoxelIndices), 'LSModelIter', params.LSModelIter, ...
        'LSSmoothnessScaling', params.LSSmoothnessScaling, 'perturbScale', params.perturbScale, ...
        'invokeLSModelPeriod', params.invokeLSModelPeriod));

    disp('start subspace move'); 
    x=subspaceMove(opts);
    x=reshape(x, size(voxelMask));
end

function xUnsharp=unsharpEnhancement(x0, voxelMask)
    xUnsharp=x0;
    H=fspecial('unsharp');
    H=H*0.1;H(2, 2)=-(sum(H(:, 1))+sum(H(:, 3))+H(1, 2)+H(3, 2))+1;
    for i=1:size(x0, 3)
        xUnsharp(:, :, i)=imfilter(x0(:, :, i), H);
    end
    xUnsharp=xUnsharp.*voxelMask;
end

function x=waveletDenoising(V, x0, sensitivities, G, data, objParams, nbrTermWeights, fullNbrMat) 
    V0coefs=LSApproxModel(V, x0, sensitivities, G, data, objParams, nbrTermWeights, fullNbrMat, objParams(2), 80, 0.6);    
    x=x0(:)+V*V0coefs;x=reshape(x, size(x0));    
end

function [segVoxelIndices, segMap]=refineSPMSegMap(segMap)
    allLabels=unique(segMap(:));allLabels(1)=[];    
    lbCount=length(allLabels);
    segVoxelIndices=cell(lbCount, 1);
    for i=1:lbCount
        segMap(segMap==allLabels(i))=i;
        segVoxelIndices{i}=find(segMap==i);
    end    
end
