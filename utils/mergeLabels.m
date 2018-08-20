function L=mergeLabels(L, labelCount)
    labels=unique(L(:));labels(1)=[];
    if length(labels)<=labelCount
        return;
    end
    components=cell(length(labels), 1);
    componentSizes=zeros(length(labels), 1);
    for i=1:length(labels)
        components{i}=find(L==labels(i));
        componentSizes(i)=length(components{i});
    end
    [~, idx]=sort(componentSizes);
    for i=1:length(componentSizes)-labelCount
        L=mergeComponent(components{idx(i)}, L);
    end
end

function L=mergeComponent(component, L)
    merged=false;
    label=L(component(1));
    direction=ceil(rand()*6);
    trialCount=1;
    while ~merged && trialCount<=6
        direction=mod(direction, 6)+1;
        curIdx=component(ceil(rand()*length(component)));
        [r, c, s]=ind2sub(size(L), curIdx);
        outBound=false;
        while ~outBound
            [r, c, s, outBound]=directionalMove(direction, r, c, s, L);
            if outBound
                break;
            end
            newLabel=L(r, c, s);
            if newLabel~=label && newLabel>=0
                L(component)=newLabel;
                merged=true;
                break;
            elseif newLabel<0
                break;
            end
        end
        trialCount=trialCount+1;
    end
    if ~merged
        L(component)=0;
    end
end

function [r, c, s, outBound]=directionalMove(direction, r, c, s, L)
    if direction==1
        r=r+1;
        outBound=r>size(L, 1);
    elseif direction==2
        c=c+1;
        outBound=c>size(L, 2);
    elseif direction==3
        s=s+1;
        outBound=s>size(L, 3);
    elseif direction==4
        r=r-1;
        outBound=r<1;
    elseif direction==5
        c=c-1;
        outBound=c<1;
    elseif direction==6
        s=s-1;
        outBound=s<1;
    end
end