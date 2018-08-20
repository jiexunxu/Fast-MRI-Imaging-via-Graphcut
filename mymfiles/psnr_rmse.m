function [psnr_noref, rmse_noref, psnr, rmse] = psnr_rmse(I, Iref)
% finds psnr and if required rmse values of different reconstructed images 
% passed parameters - I or filename, refimg or filename
% if no params passed, a GUI is presented to choose files

medfiltwin = [];  %21;% window for median filtering to correct bias - used to compute rmse; [] to disable local bias corr
roi_flg = 0;      % 1 enables ROI selection
del = 3;              % removes this much from image margins when computing noise etc
thisdir = pwd;
if nargin==0
    % datadir = pwd;
    % if ~isempty(datadir),  cd(datadir);  end
    [user_datafile, user_datapath] = uigetfiles('*.*', 'Pick input image data file, and if rmse required, also the reference file in order (I_flname, ref_flname)');
    if (isequal(user_datafile,0) | isequal(user_datapath,0))
        disp('no file(s) selected');
        cd(thisdir);
        return
    else
        data_name = user_datafile;
        datadir = user_datapath;   
        cd(datadir);
        if iscell(user_datafile)
            listlen = length(user_datafile);
            if listlen ~= 2
                disp('select at most 2 image files');  
                cd(thisdir);
                return
            end
            % wrap list by one - needed for some reason...
           tmp = user_datafile;
           for kk = 2:listlen, data_name{kk-1} = tmp{kk};  end
           data_name{listlen} = tmp{1};
           clear tmp;
           disp(data_name(:));
           I_flname = data_name{1};
           Iref_flname = data_name{2};       
           Iref = (read_img(Iref_flname)); 
           I = (read_img(I_flname)); 
        else
            I_flname = data_name;
            I = (read_img(I_flname)); 
        end        
    end
elseif nargin==1
    if ischar(I)
        I_flname = I;
        I = (read_img(I_flname)); 
    end        
elseif nargin==2
    if ischar(I)
        I_flname = I;
        I = (read_img(I_flname)); 
    end        
    if ischar(Iref)
        Iref_flname = I;
        Iref = (read_img(Iref_flname)); 
    end        
end
cd(thisdir);

if ~exist('Iref', 'var')
    I = I*255/max(abs(I(:)));
    % compute psnr from I only
    [psnr_noref, rmse_noref] = get_psnr_noref(I);
else
    sref = size(Iref); sI = size(I);
    % make same size
    minsiz = min([sref; sI]);
    if ndims(I) ~= ndims(Iref)
        disp(' 2 arrays of different dimensions!');
    end
    if ~isequal(sref, sI), disp('warning: 2 arrays of different size - snr result might be erroneous!');  end
    if ndims(I)==2
        Iref = Iref(end - minsiz(1)+1:end, end - minsiz(2)+1:end);
        I = I(end - minsiz(1)+1:end, end - minsiz(2)+1:end);
    elseif ndims(I)==3
        if ~isequal(sref(3), sI(3)), disp('Error: 2 3D arrays whose number of slices is unequal - cannot find snr!');  end
        Iref = Iref(end - minsiz(1)+1:end, end - minsiz(2)+1:end, end - minsiz(3)+1:end);
        I = I(end - minsiz(1)+1:end, end - minsiz(2)+1:end, end - minsiz(3)+1:end);
    elseif ndims(I)==4
        if ~isequal(sref(4), sI(4)), disp('Error: 2 4D arrays whose number of slices or phases is unequal - cannot find snr!');  end
        Iref = Iref(end - minsiz(1)+1:end, end - minsiz(2)+1:end, end - minsiz(3)+1:end, end - minsiz(4)+1:end);
        I = I(end - minsiz(1)+1:end, end - minsiz(2)+1:end, end - minsiz(3)+1:end, end - minsiz(4)+1:end);
    end        
    % equalize max
    mxref = max(abs(Iref(:))); mxI = max(abs(I(:)));
    %     figure;  imagesc([Iref, I]); colormap(gray)
    Iref = Iref*255/mxref;
    I = I*255/mxI;
    % compute psnr from I only
    [psnr_noref, rmse_noref] = get_psnr_noref(I);
    % compute psnr and rmse from I-Iref
    [psnr, rmse] = get_psnr_rmse(I, Iref);
end    
    
%     if roi_flg
%         figure;
%         [BW, xi, yi] = roipoly(I0/max(max(I0)));
%         %figure, imshow(I0), figure, imshow(BW)
%         rectangle('Position', [xi(1), yi(1), xi(2) - xi(1), yi(2) - yi(1)])
%         d0 = d0(yi(1):yi(2), xi(1):xi(2));
%         dgc = dgc(yi(1):yi(2), xi(1):xi(2));
%         rectcoords = [yi(1:2), xi(1:2)],
%     end
    

    
%% Internal functions

function [psnr, rmse] = get_psnr_rmse(x, x0)
    x0 = single(x0);
    x = single(x); 
    k = 0;
    for shifti = -1:1
        for shiftj = -1:1
            k = k+1;
            x1 = shift(x, shifti, shiftj);
            d = x1-x0;
            if isempty(medfiltwin)
                d = d(del+1:end-del,del+1:end-del);
                xorig = x0(del+1:end-del,del+1:end-del);
                dmean = mean(d(:));
            else
                if ndims(d)==2
                    d = d(del+1:end-del,del+1:end-del);
                    xorig = x0(del+1:end-del,del+1:end-del);
                    dmean = medfilt2(d, [medfiltwin, medfiltwin]);
                else
                    szx = size(d);
                    dd = reshape(d, [szx(1), prod(szx(2:end))]);
                    dmean = medfilt2(dd, [medfiltwin, medfiltwin]);
                    clear dd;
                    dmean = reshape(dmean, szx);
                    d = d(del+1:end-del,del+1:end-del,:);
                    dmean = dmean(del+1:end-del,del+1:end-del,:);
                    xorig = x0(del+1:end-del,del+1:end-del,:);
                end
            end
            d = d - dmean;
            rmsevec(k) = (norm(d(:))/ sqrt(numel(d))) / (norm(xorig(:))/ sqrt(numel(xorig)));    
            psnrvec(k) = -20*log10(rmsevec(k) / max(abs(xorig(:))));      
        end
    end
    psnr = max(psnrvec);
    rmse = min(rmsevec);
end   % of get_psnr_rmse

function img = shift(im, i, j)
    if ndims(im)==2
        img = im;
        [m,n] = size(img);
        if i > 0
            img = [zeros(i,n); img(1:end-i, :)];
        elseif i < 0
            img = [img(abs(i)+1:end, :); zeros(abs(i), n)];
        end;
        [m,n] = size(img);
        if j > 0
            img = [zeros(m,j), img(:,1:end-j)];
        elseif j < 0
            img = [img(:, abs(j)+1:end), zeros(m, abs(j))];
        end;
    elseif ndims(im)>2
        sz = size(im);
        img = reshape(im, [sz(1), sz(2), prod(sz(3:end))]);
        [m,n, p] = size(img);
        if i > 0
            img = cat(1, zeros(i,n,p), img(1:end-i, :,:));
        elseif i < 0
            img = cat(1, img(abs(i)+1:end, :,:), zeros(abs(i), n,p));
        end;
        [m,n,p] = size(img);
        if j > 0
            img = cat(2, zeros(m,j,p), img(:,1:end-j,:));
        elseif j < 0
            img = cat(2, img(:, abs(j)+1:end,:), zeros(m, abs(j),p));
        end;
        img = reshape(img, sz);
    end
end   % of shift

function [psnr, rmse] = get_psnr_noref(x)
% finds psnr without reference image, by extracting noise from image
        % get hpf(x)
        szx = size(x);
        if ndims(x)>2
            x = reshape(x, [szx(1), szx(2), prod(szx(3:end))]);
        end
        [nrows, ncols, nslices] = size(x);
        sigedge = 1.5;
        sighp = 2;
        h = fspecial('Gaussian', 5, sighp);
        h = h/sum(sum(h));
        h = max(max(h)) - h;

        % 1) edge det
        for i=1:nslices
            q = edge(abs(x(:,:,i)), 'canny', [], sigedge);
            edg(:,:,i) = bwmorph(q, 'dilate');
        end

        %2) hig pass filtering
        for i=1:nslices
            d(:,:,i) = conv2(x(:,:,i), h, 'same');
        end

        %3) make sure edges and image margins are deleted from HPF
        d(edg) = 0;   
        d = d(del+1:end-del,del+1:end-del,:);
        xorig = x(del+1:end-del,del+1:end-del,:);
        mx = max(abs(xorig(:)));
        rmse = (norm(d(:))/ sqrt(numel(d))) / (norm(xorig(:))/ sqrt(numel(xorig)));    
        psnr = -20*log10(rmse/mx);      
end

end % of main