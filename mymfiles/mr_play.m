%function mr_play;
% reads a MR angio series (in k-space), transforms into img spce
% shows the series as a movie

% % data = 'Eberheart_Annie';
% data = 'Bruckner_Donald';
% data = 'Cupid_Charmaine';
% load(data);

%dataformat = 'cell';
dataformat = 'array';

if strcmp(dataformat, 'cell')
    
    for i = 1:length(Series)
        %C{i} = ifft2(doXchop(Series{i}));
        C{i} = ifft2(doYchop(Series{i}));
    end

    % % displays movie of the angio frames in image-space
    % figure;
    % for i = 1:length(C)
    % imshow(abs(C{i}));
    % M(i) = getframe;
    % end;
    % 
    % movie(M, 2, 2);
    % 
    % displays movie of the angio difference in image-space
    figure;
    for i = 2:length(C)
        q = C{i} - C{i-1};
    imagesc(abs(q));  colormap(gray);  title(sprintf('frame %d', i));
    dummy = input('dummy');
    M(i) = getframe;
    end;
    movie(M(2:end), 2, 2);
elseif strcmp(dataformat, 'array')
    [m,n,c] = size(C);
    for i = 1:c
        %C(:,:,i) = ifft2(doXchop(C(:,:,i)));
        C(:,:,i) = ifft2(doYchop(C(:,:,i)));
    end

    % % displays movie of the angio frames in image-space
    % figure;
    % for i = 1:c
    % imshow(abs(C(:,:,i)));
    % M(i) = getframe;
    % end;
    % 
    % movie(M, 2, 2);
    % 
    % displays movie of the angio difference in image-space
    figure;
    for i = 2:c
        q = C(:,:,i) - C(:,:,i-1);
    imagesc(abs(q));  colormap(gray);  title(sprintf('frame %d', i));
    dummy = input('dummy');
    M(i) = getframe;
    end;
    movie(M(2:end), 2, 2);
end

% % Undersampling/aliasing example
% R = 3;  % undersampling factor
% figure;
% for i = 1:length(Series)
%   q = fft2(C{i});
%   q = q(1:R:end, :);   % undersampling in k_y direction
%   CC{i} = ifft2(q);
%   imshow(abs(CC{i}));
%   MM(i) = getframe;
% end
% 
% movie(MM);
  






