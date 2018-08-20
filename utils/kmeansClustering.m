function [idx, L]=kmeansClustering(x, k)
    dataMat=zeros(numel(x), 4);
    [I1, I2, I3]=ind2sub(size(x), numel(x));
    dataMat(:, 1)=I1/max(I1);dataMat(:, 2)=I2/max(I2);
    dataMat(:, 3)=I3/max(I3);dataMat(:, 4)=0*x(:)/max(x(:));
    idx=kmeans(dataMat, k, 'emptyaction', 'drop', 'start', 'cluster');
    idx=reshape(idx, size(x));
    L=label2rgb(idx);
end