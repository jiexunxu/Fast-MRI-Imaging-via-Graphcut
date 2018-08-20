function Vout=sparsifyDenseSubspace(V, segVoxelIndices, voxelMask)    
    Vrows=zeros(nnz(voxelMask)*size(V, 2), 1);
    Vcols=zeros(length(Vrows), 1);Vvals=zeros(length(Vrows), 1);
    Vcounter=1;VidxCounter=1;
    for i=1:length(segVoxelIndices)
        segVoxelIdx=segVoxelIndices{i};
        if ~isempty(segVoxelIdx)
            for j=1:size(V, 2)                                     
                Vidx=VidxCounter:VidxCounter+length(segVoxelIdx)-1;
                Vrows(Vidx)=segVoxelIdx;Vcols(Vidx)=Vcounter;Vvals(Vidx)=V(segVoxelIdx, j);
                Vcounter=Vcounter+1;
                VidxCounter=VidxCounter+length(segVoxelIdx);
            end
        end
    end
    Vrows(VidxCounter:length(Vrows))=[];
    Vcols(VidxCounter:length(Vcols))=[];
    Vvals(VidxCounter:length(Vvals))=[];
    Vout=sparse(Vrows, Vcols, Vvals, numel(voxelMask), Vcounter-1);  
end