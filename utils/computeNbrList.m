function [nbrList, nbrAdjList, nbrWeights]=computeNbrList(voxelMask, varargin)   
    nbrCount=6;
    [numRows, numCols, numSlices]=size(voxelMask);
    edgeMap=zeros(numRows, numCols, numSlices);
    outputNbrWeights=false;
    if nargin==2
        edgeMap=varargin{1};
    elseif nargin==3
        % edgeWeights(1) is edge - nonedge voxel pair weights, edgeWeights(2) is
        % edge - edge voxel pair weights
        edgeMap=varargin{1};
        edgeWeights=varargin{2};
        outputNbrWeights=true;        
    end
    voxelMaskList=find(voxelMask(:));
    nbrList=zeros(nbrCount*length(voxelMaskList)/2, 2);
    nbrAdjList=zeros(numel(voxelMask), nbrCount);
    nbrAdjListCounters=ones(numel(voxelMask), 1);
    nbrListCounter=1;
    if outputNbrWeights
        nbrWeights=zeros(nbrCount*length(voxelMaskList)/2, 1);
    end
    for i=1:length(voxelMaskList)
        vid=voxelMaskList(i);
        [r, c, s]=idx2rcs(vid);
        if nbrCount>=6 && (~edgeMap(r, c, s) || outputNbrWeights)
            if s<numSlices && voxelMask(r, c, s+1) && (~edgeMap(r, c, s+1) || outputNbrWeights)                    
                nbrList(nbrListCounter, 1)=vid;
                nbrList(nbrListCounter, 2)=rcs2idx(r, c, s+1);
                if outputNbrWeights                   
                    if edgeMap(r, c, s) && edgeMap(r, c, s+1)
                        nbrWeights(nbrListCounter)=edgeWeights(2);
                    elseif edgeMap(r, c, s) && ~edgeMap(r, c, s+1)
                        nbrWeights(nbrListCounter)=edgeWeights(1);
                    else
                        nbrWeights(nbrListCounter)=1;
                    end
                end
                nbrListCounter=nbrListCounter+1;
                nbrAdjList(vid, nbrAdjListCounters(vid))=rcs2idx(r, c, s+1);
                nbrAdjListCounters(vid)=nbrAdjListCounters(vid)+1;               
            end
            if s>1 && voxelMask(r, c, s-1) && (~edgeMap(r, c, s-1) || outputNbrWeights)
                nbrAdjList(vid, nbrAdjListCounters(vid))=rcs2idx(r, c, s-1);
                nbrAdjListCounters(vid)=nbrAdjListCounters(vid)+1;
            end
            if c<numCols && voxelMask(r, c+1, s) && (~edgeMap(r, c+1, s) || outputNbrWeights)
                nbrList(nbrListCounter, 1)=vid;
                nbrList(nbrListCounter, 2)=rcs2idx(r, c+1, s);
                if outputNbrWeights
                    if edgeMap(r, c, s) && edgeMap(r, c+1, s)
                        nbrWeights(nbrListCounter)=edgeWeights(2);
                    elseif edgeMap(r, c, s) && ~edgeMap(r, c+1, s)
                        nbrWeights(nbrListCounter)=edgeWeights(1);
                    else
                        nbrWeights(nbrListCounter)=1;
                    end
                end
                nbrListCounter=nbrListCounter+1;
                nbrAdjList(vid, nbrAdjListCounters(vid))=rcs2idx(r, c+1, s);
                nbrAdjListCounters(vid)=nbrAdjListCounters(vid)+1;  
            end
            if c>1 && voxelMask(r, c-1, s) && (~edgeMap(r, c-1, s) || outputNbrWeights)
                nbrAdjList(vid, nbrAdjListCounters(vid))=rcs2idx(r, c-1, s);
                nbrAdjListCounters(vid)=nbrAdjListCounters(vid)+1;  
            end
            if r<numRows && voxelMask(r+1, c, s) && (~edgeMap(r+1, c, s) || outputNbrWeights)
                nbrList(nbrListCounter, 1)=vid;
                nbrList(nbrListCounter, 2)=rcs2idx(r+1, c, s);
                if outputNbrWeights
                    if edgeMap(r, c, s) && edgeMap(r+1, c, s)
                        nbrWeights(nbrListCounter)=edgeWeights(2);
                    elseif edgeMap(r, c, s) && ~edgeMap(r+1, c, s)
                        nbrWeights(nbrListCounter)=edgeWeights(1);
                    else
                        nbrWeights(nbrListCounter)=1;
                    end
                end
                nbrListCounter=nbrListCounter+1;
                nbrAdjList(vid, nbrAdjListCounters(vid))=rcs2idx(r+1, c, s);
                nbrAdjListCounters(vid)=nbrAdjListCounters(vid)+1;  
            end     
            if r>1 && voxelMask(r-1, c, s) && (~edgeMap(r-1, c, s) || outputNbrWeights)
                nbrAdjList(vid, nbrAdjListCounters(vid))=rcs2idx(r-1, c, s);
                nbrAdjListCounters(vid)=nbrAdjListCounters(vid)+1;  
            end
        end
        
        if nbrCount>=18 && (~edgeMap(r, c, s) || outputNbrWeights)
            if c<numCols && s<numSlices && voxelMask(r, c+1, s+1)  && ~edgeMap(r, c+1, s+1)
                nbrList(nbrListCounter, 1)=vid;
                nbrList(nbrListCounter, 2)=rcs2idx(r, c+1, s+1);
                if outputNbrWeights                   
                    if edgeMap(r, c, s) && edgeMap(r, c+1, s+1)
                        nbrWeights(nbrListCounter)=edgeWeights(2);
                    elseif edgeMap(r, c, s) && ~edgeMap(r, c+1, s+1)
                        nbrWeights(nbrListCounter)=edgeWeights(1);
                    else
                        nbrWeights(nbrListCounter)=1;
                    end
                end
                nbrListCounter=nbrListCounter+1;
            end
            if r<numRows && c<numCols && voxelMask(r+1, c+1, s) && ~edgeMap(r+1, c+1, s)
                nbrList(nbrListCounter, 1)=vid;
                nbrList(nbrListCounter, 2)=rcs2idx(r+1, c+1, s);
                if outputNbrWeights                   
                    if edgeMap(r, c, s) && edgeMap(r+1, c+1, s)
                        nbrWeights(nbrListCounter)=edgeWeights(2);
                    elseif edgeMap(r, c, s) && ~edgeMap(r+1, c+1, s)
                        nbrWeights(nbrListCounter)=edgeWeights(1);
                    else
                        nbrWeights(nbrListCounter)=1;
                    end
                end
                nbrListCounter=nbrListCounter+1;
            end
            if r<numRows && s<numSlices && voxelMask(r+1, c, s+1) && ~edgeMap(r+1, c, s+1)
                nbrList(nbrListCounter, 1)=vid;
                nbrList(nbrListCounter, 2)=rcs2idx(r+1, c, s+1);
                if outputNbrWeights                   
                    if edgeMap(r, c, s) && edgeMap(r+1, c, s+1)
                        nbrWeights(nbrListCounter)=edgeWeights(2);
                    elseif edgeMap(r, c, s) && ~edgeMap(r+1, c, s+1)
                        nbrWeights(nbrListCounter)=edgeWeights(1);
                    else
                        nbrWeights(nbrListCounter)=1;
                    end
                end
                nbrListCounter=nbrListCounter+1;
            end
        end
        
        if nbrCount>=26
            if r<numRows && s<numSlices && c<numCols && voxelMask(r+1, c+1, s+1) && ~edgeMap(r+1, c+1, s+1)
                nbrList(nbrListCounter, 1)=vid;
                nbrList(nbrListCounter, 2)=rcs2idx(r+1, c+1, s+1);
                if outputNbrWeights                   
                    if edgeMap(r, c, s) && edgeMap(r+1, c+1, s+1)
                        nbrWeights(nbrListCounter)=edgeWeights(2);
                    elseif edgeMap(r, c, s) && ~edgeMap(r+1, c+1, s+1)
                        nbrWeights(nbrListCounter)=edgeWeights(1);
                    else
                        nbrWeights(nbrListCounter)=1;
                    end
                end
                nbrListCounter=nbrListCounter+1;
            end
        end
        
    end
    nbrList(nbrListCounter:size(nbrList, 1), :)=[];
    if outputNbrWeights
        nbrWeights(nbrListCounter:size(nbrWeights, 1))=[];
    end
    
    function [r, c, s]=idx2rcs(idx)
        s=ceil(idx/(numRows*numCols));
        idx=idx-(s-1)*numRows*numCols;
        c=ceil(idx/numRows);
        r=idx-(c-1)*numRows;
    end

    function idx=rcs2idx(r, c, s)
        idx=r+(c-1)*numRows+(s-1)*numRows*numCols;
    end
end