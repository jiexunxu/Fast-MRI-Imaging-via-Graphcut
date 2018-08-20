function [read_success, found_dicom] = preview_dicoms()
user_dir = pwd;
D = dir(pwd);
imgind = find(~[D.isdir]);
found_dicom = false;
read_success = false;
try
    q = dicomread(D(imgind(1)).name);
    disp('successfully read DICOMS from this directory');
    read_success = true;
catch
    disp('DICOM read not successful - I will now show previews of subfolders');
    dirind = find([D.isdir]);
    ndirs = length(dirind); 
    if ndirs < 3, disp('no subdirectories found! - cant continue'); end
    k=1;
    fig = figure;
    for i = 1:length(dirind)
        if strcmp(D(dirind(i)).name, '.') | strcmp(D(dirind(i)).name, '..')
            continue;
        end
        D(dirind(i)).name,
        cd(D(dirind(i)).name);
        DD = dir(pwd),
        imgind2 = find(~[DD.isdir]);
        nimages = length(imgind2);
        try
            q = dicomread(DD(imgind2(1)).name);
            try 
                qq = dicominfo(DD(imgind2(1)).name);
                seriesdesc = qq.SeriesDescription;
            catch
                seriesdesc = 'not found';
            end
            if isempty(q), disp('the first file in directory is not a readable dicom'); continue; end
            thisdir = pwd;
            if length(thisdir) > 25
                thisdir = thisdir(end-25:end);
            end
            subplot(ceil(sqrt(ndirs-2)), ceil(sqrt(ndirs-2)), k); imagesc(q); title(sprintf('dir: %s, #images: %d, Series: %s', thisdir, nimages, seriesdesc))
            k = k+1;
            found_dicom = true;
            cd(user_dir);
        catch
            cd(user_dir);
            continue;
        end
    end
    figure(fig);
end
cd(user_dir);
disp('Press something...');
pause;
