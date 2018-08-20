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