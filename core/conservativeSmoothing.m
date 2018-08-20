function img=conservativeSmoothing(img, w)
    img=abs(img);
    W=size(img, 1);H=size(img, 2);
    for i=w+1:W-w
        for j=w+1:H-w
            nbr=img(i-w:i+w, j-w:j+w);nbr=nbr(:);
            nbr((length(nbr)+1)/2)=[];
            if img(i, j)<min(nbr)
                img(i, j)=min(nbr);
            elseif img(i, j)>max(nbr)
                img(i, j)=max(nbr);
            end
        end
    end
end