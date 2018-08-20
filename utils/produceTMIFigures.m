function produceTMIFigures(ref, SENSE, CS, Graphcut)
    x0=abs(ref);x1=abs(SENSE);x2=abs(CS);x3=abs(Graphcut);
    x0=x0/median(x0(:));x1=x1/median(x1(:));
    x2=x2/median(x2(:));x3=x3/median(x3(:));
    x0=resizem(x0, [512 512], 'bicubic');
    x1=resizem(x1, [512 512], 'bicubic');
    x2=resizem(x2, [512 512], 'bicubic');
    x3=resizem(x3, [512 512], 'bicubic');
    visualizeImages([2, 2], x1, x2, x3, x0);
    dif1=abs(x0-x1);dif2=abs(x0-x2);dif3=abs(x0-x3);
    m=max([dif1(:); dif2(:); dif3(:)]);
    dif1=uint8(dif1*255/m);dif2=uint8(dif2*255/m);dif3=uint8(dif3*255/m);
    figure;
    subaxis(2, 2, 1);colormap(gray);imagesc(dif1);axis off;
    subaxis(2, 2, 2);colormap(gray);imagesc(dif2);axis off;
    subaxis(2, 2, 3);colormap(gray);imagesc(dif3);axis off;    
end