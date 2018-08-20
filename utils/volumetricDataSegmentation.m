function L=volumetricDataSegmentation(x, edgeMap)  
    x=gaussianFilter3D(x, [7 7], 7);   
    x=preprocessingData(x, edgeMap);
    L=watershed(x, 6);
  %  L=int32(L(2:m+1, 2:n+1, 2:p+1));
    L(x==0)=-1;
%    L=mergeLabels(L, labelCount);
end

function x=gaussianFilter3D(x, w, std)
    for i=1:size(x, 3)
        x(:, :, i)=imfilter(x(:, :, i), fspecial('gaussian', [w(1), w(2)], std));
    end
    %{
    [Fx, Fy, Fz]=meshgrid(-w(1):w(1), -w(2):w(2), -w(3):w(3));
    arg=-(Fx.*Fx+Fy.*Fy+Fz.*Fz)/(2*std*std);
    gaussian=exp(arg);
    gaussian=gaussian/sum(gaussian(:));
    x=imfilter(x, gaussian, 0);
    %}
end

function gradmag=preprocessingData(x, edgeMap)
%{
    sobel3dBase=[1 2 1; 2 4 2; 1 2 1];
    sobelX=zeros(3, 3, 3);sobelY=zeros(3, 3, 3);sobelZ=zeros(3, 3, 3);
    sobelX(1, :, :)=-sobel3dBase;sobelX(3, :, :)=sobel3dBase;
    sobelY(:, 1, :)=-sobel3dBase;sobelY(:, 3, :)=sobel3dBase;
    sobelZ(:, :, 1)=-sobel3dBase;sobelZ(:, :, 3)=sobel3dBase;
    Ix=imfilter(x, sobelX, 0);Iy=imfilter(x, sobelY, 0);Iz=imfilter(x, sobelZ, 0);
%}
    [m, n, p]=size(x);
    Ix=zeros(m, n, p);Iy=zeros(m, n, p);
    for i=1:size(x, 3)
        Ix(:, :, i)=imfilter(x(:, :, i), fspecial('sobel'));
        Iy(:, :, i)=imfilter(x(:, :, i), fspecial('sobel')');
    end
    gradmag=sqrt(Ix.^2+Iy.^2);
    gradmag(edgeMap==1)=gradmag(edgeMap==1)*1.3;
    gradmag=gradmag-min(gradmag(:));gradmag=gradmag/max(gradmag(:));
    stepValues=linspace(0, 1, 100+1);
    for i=1:length(stepValues)-1
        lower=stepValues(i);upper=stepValues(i+1);
        gradmag(gradmag>=lower & gradmag<=upper)=lower;
    end
end