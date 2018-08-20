function [v_idx, vecCount]=pickVecs(V, vecCount, totalSupport, iter, nullVecSelectionMethod)
    N=size(V, 2);
    if strcmpi(nullVecSelectionMethod, 'sequential')
        start_idx=mod((iter-1)*vecCount, size(V, 2))+1;        
        end_idx=mod(iter*vecCount-1, size(V, 2))+1;        
        if start_idx>=end_idx
            v_idx=[start_idx:size(V, 2) 1:end_idx];
        else
            v_idx=start_idx:end_idx;
        end        
    elseif strcmpi(nullVecSelectionMethod, 'random')
        if vecCount<N/2
            vecs=zeros(N, 1);
            k=0;
            while k<vecCount
                nextIdx=ceil(rand*N);
                if vecs(nextIdx)==0
                    vecs(nextIdx)=1;
                    k=k+1;
                end
            end
        else
            vecs=ones(N, 1);
            k=N;
            while k>vecCount
                nextIdx=ceil(rand*N);
                if vecs(nextIdx)==1
                    vecs(nextIdx)=0;
                    k=k-1;
                end
            end
        end
        v_idx=find(vecs); 
    elseif strcmpi(nullVecSelectionMethod, 'graphcut')
        %{
        attempts=0;
        supportedList=false(size(V, 1), 1);
        v_idx=zeros(vecCount, 1);
        v_idx_counter=1;
        while attempts<vecCount*2 && v_idx_counter<=vecCount
            attempts=attempts+1;
            nextVec=min(abs(ceil(randn*size(V, 2)/3))+1, size(V, 2));
            v=V(:, nextVec);
            numSupports=sum(supportedList & (v~=0));
            if numSupports==0
                supportedList(v~=0)=true;
                v_idx(v_idx_counter)=nextVec;
                v_idx_counter=v_idx_counter+1;
            else
                if rand<1/10^numSupports
                    v_idx(v_idx_counter)=nextVec;
                    v_idx_counter=v_idx_counter+1;
                    supportedList(v~=0)=true;
                end
            end
        end        
        v_idx(v_idx_counter:length(v_idx))=[];
        %}
        
        v_idx=zeros(vecCount, 1);
        v_idx_counter=1;
        supportedList=zeros(size(V, 1), 1);
        shuffledList=randperm(size(V, 2));
        for i=1:size(V, 2)
            nextVec=shuffledList(i);
            v=V(:, nextVec);
            numSupports=sum(supportedList & (v~=0));
            
            if numSupports==0 || rand<=0/10^numSupports
         %   if max(supportedList-10*ones(length(supportedList), 1))<=0
                supportedList(v~=0)=true;
                supportedList(v~=0)=supportedList(v~=0)+1;
                v_idx(v_idx_counter)=nextVec;
                v_idx_counter=v_idx_counter+1;
            end
        end
        v_idx(v_idx_counter:length(v_idx))=[];
        
        fprintf('%d (%d percent) voxels are supported\n', sum(supportedList), sum(supportedList)*100/totalSupport);
    end
end