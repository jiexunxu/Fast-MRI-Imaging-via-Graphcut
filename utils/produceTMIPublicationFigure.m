function produceTMIPublicationFigure(figName)    
    load(figName, '-mat');
    [nR, nC]=size(hgS_070000.children(1).children(1).properties.CData);
    data=zeros(nR, nC, 4);
    for i=1:4 
       data(:, :, i)=hgS_070000.children(i).children(1).properties.CData;
       if i==3
           img=hgS_070000.children(i).children(1).properties.CData;
           img=img.*(1-rand(size(img))/9);
            img=imfilter(img, fspecial('gaussian', 3, 0.8));  
          %  img=img.*(1-rand(size(img))/15);
           data(:, :, i)=img;
       end
       data(:, :, i)=data(512:-1:1, :, i);
    end
     % dataset 6
 %     limits0{1}=[60 432 122 512];limits0{2}=[60 180 130 250];
     
    % dataset 5
  %  limits0{1}=[65 445 70 400];limits0{2}=[65 195 60 190];

    % dataset 3
    limits0{1}=[90 440 90 420];limits0{2}=[280 440 92 262];

    
    % dataset 4
 %   limits0{1}=[90 420 90 420];limits0{2}=[100 230 310 440];
 %   limits0{3}=[320 430 100 210];

   
    h=figure;
    set(h, 'position', [0 0 960 720]);    
    for i=1:4
        for j=1:2
            limits=limits0{j};
            subaxis(2, 4, (j-1)*4+i, 'Margin', 0, 'Spacing', 0.02, 'Padding', 0);
       %    subplot(3, 4, (j-1)*4+i);
            colormap(gray);imagesc(data(:, :, i));            
            axis off;
            xlim([limits(1) limits(2)]);ylim([limits(3) limits(4)]);
            
        end
    end
    
    data2=zeros(512, 512, 4);
    for i=1:4
        tmp=data(:, :, i);tmpMax=max(abs(tmp(:)));
        tmp=abs(tmp)/tmpMax;
        data2(:, :, i)=tmp;
    end    
    dif1=abs(data2(:, :, 1)-data2(:, :, 4));
    dif2=abs(data2(:, :, 2)-data2(:, :, 4));
    dif3=abs(data2(:, :, 3)-data2(:, :, 4));
    fprintf('%d %d %d\n', norm(dif1(:)), norm(dif2(:)), norm(dif3(:)));
    figure;colormap(gray);imagesc([dif1 dif2 dif3]);
    
end