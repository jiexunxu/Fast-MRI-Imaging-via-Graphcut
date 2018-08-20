function params=subspaceGraphcutParameters
    platform=computer;
    if strcmpi(platform, 'PCWIN64')
        params=windowsLaptop_params();
    elseif strcmpi(platform, 'GLNXA64')
        params=p41_params();
    end
end

function params=windowsLaptop_params
    params.initializationMethod='SENSE';
    
    params.SENSEweight=0.5;
    
    params.L1MaxIter=10;
    params.L1Sigma=1e-9;
    params.L1SmoothnessWeight=0;
           
    params.dataTermWeight=1;
    params.EdEsRatio=0.5;
    params.truncationFactor=0.1;
    params.smoothnessPower=0.8;
        
    params.BFsimgas{1}=linspace(0.6, 1.4, 10);
    params.BFsigmaWeight=0.011;
    params.use3DBF=false;
    
    params.waveletDecompLevel=4;
    params.waveletFilter='db3';
    params.combineSameLvlWavelet=false;
    
    params.segmentationMethod='notSPM8Volumetrics';
    params.graphSegmentationMinSizeRatio=0.0015;
    params.whiteMatterSegmentationMinSizeRatio=0.03;
    
    params.edgeMapThreshold=0.01;
    
    params.graphcutIterCount=20;
        
    params.LSModelIter=40;
    params.LSSmoothnessScaling=0.9;
    params.perturbScale=0;
    params.invokeLSModelPeriod=3;
end

function params=p41_params
    params.initializationMethod='CS';
    
    params.SENSEweight=0.3;
    
    params.L1MaxIter=80;
    params.L1Sigma=1e-9;
    params.L1SmoothnessWeight=0.003;
           
    params.dataTermWeight=1;
    params.EdEsRatio=0.3;
    params.truncationFactor=0.3;
    params.smoothnessPower=0.8;
        
    params.BFsimgas{1}=linspace(2, 4, 5);
    params.BFsigmaWeight=0.011;
    params.use3DBF=true;
    
    params.segmentationMethod='SPM8Volumetrics';
    params.graphSegmentationMinSizeRatio=7e-4;
    params.whiteMatterSegmentationMinSizeRatio=0.03;
    
    params.edgeMapThreshold=0.01;
    
    params.graphcutIterCount=30;
        
    params.LSModelIter=50;
    params.LSSmoothnessScaling=0.5;
    params.perturbScale=0.05;
    params.invokeLSModelPeriod=3;
end