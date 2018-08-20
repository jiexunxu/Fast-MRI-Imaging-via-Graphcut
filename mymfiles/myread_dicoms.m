function [dcm_vol, volsize, voxsize, ndicoms, info] = myread_dicoms(data_dirname)
% reads in all dicoms from a folder, 
% and rearranges them into a 3d or 4d volume according to 
% the table position and acquisition time

thisdir = pwd;
if nargin==0
    data_dirname = thisdir;
end
info = [];
% ask for input directory 
user_path = uigetdir(data_dirname, 'Select directory containing dicoms (CANCEL to accept suggested)');
if ~isequal(user_path,0)    
    data_dirname = user_path;
end
cd(data_dirname);
[read_success, found_dicom] = preview_dicoms;
if ~read_success
    if found_dicom
        % ask again for input directory 
        user_path = uigetdir(data_dirname, 'Try again: Select directory containing dicoms of CE_MRI (CANCEL to accept suggested)');
        if ~isequal(user_path,0)    
            data_dirname = user_path;
        end
        cd(data_dirname);
    else
        error('selectedt folder does not contain subfolders with readable DICOMS!');
    end
end

D = dir(pwd);
imgind = find(~[D.isdir]);

for i = 1:length(imgind)
    [val, err] = dcm_get_property(D(imgind(i)).name, {'InstanceNumber', 'AcquisitionTime', 'Filename', 'SliceThickness', 'PixelSpacing', 'ImagePositionPatient'});
    if err, try info = dicominfo(D(imgind(i)).name); end; end
        if isempty(val{1}), 
            try
                ins_all(i) = info.InstanceNumber;
            catch
                ins_all(i) = i; 
            end
        else
            ins_all(i) = str2num(val{1}); 
        end
        if isempty(val{2}), 
            try
                acq_all(i) = info.AcquisitionTime;
            catch
                acq_all(i) = 1; 
            end
        else
            acq_all(i) = str2num(val{2}); 
        end
        if isempty(val{3}), 
            try Filename_all{i} = info.Filename;
            catch Filename_all{i} = D(imgind(i)).name; end
        else
            Filename_all{i} = val{3}; 
        end
        if ~isempty(val{6})
            [ss1, dd1] = strtok(val{6}, '\d'); [ss1, dd1] = strtok(dd1(2:end), '\d'); 
            zpos(i) = str2num(dd1(2:end));
        else
            try tt = info.ImagePositionPatient; zpos(i) = tt(3); end
        end
        is_dcm(i) = isdicom(D(imgind(i)).name);
end
if ~isempty(val{5})
    [ss, dd] = strtok(val{5}, '\d');
    xres = str2num(ss); 
    [ss, dd] = strtok(dd(2:end), '\d');
    yres = str2num(ss);
    if ~isempty(dd), 
        zres =  str2num(dd(2:end)); 
    else
        szpos = sort(zpos);
        del = abs(szpos(2:end) - szpos(1:end-1));
        del = del(del>0);
        zres = min(del);
    end
else
    try
        tt = info.PixelSpacing; 
        xres = tt(1); yres = tt(2); 
        if length(tt)==3
            zres = tt(3);
        else
            szpos = sort(zpos);
            del = abs(szpos(2:end) - szpos(1:end-1));
            del = del(del>0);
            zres = min(del);
        end
    catch
        disp('voxel resolution could not be reliably determined - check header');
    end
end
voxsize = [xres, yres, zres];
ndicoms = length(find(is_dcm));

zunique = sort(unique(zpos));
for i = 1:length(zunique)
    ind = find(zpos == zunique(i));
    [y, ind2] = sort(acq_all(ind));
    %[ind, 1e-5*y, ind(ind2)],
    indextable(i,:) = ind(ind2);
end

nslices = size(indextable, 1);
nacq = size(indextable, 2);
if ~(length(imgind) == nslices*nacq)
    disp('WARNING: discrepancy between noof files, nslices and nacquisitions, by zpostion and acqtime');
    nslices = floor(length(imgind)/nacq);
    nooffiles = length(imgind), mxInstanceNumber = max(ins_all), nacq, nslices,
end
q = single(dicomread(Filename_all{1}));
[nrows, ncols] = size(q);
dcm_vol = zeros(nrows, ncols, nslices, nacq, 'single');

for i = 1:nslices
    for j = 1:nacq
        try
            q = single(dicomread(Filename_all{indextable(i,j)}));
            dcm_vol(:,:,i,j) = q;
        catch
            % dont read this file - not dicom or not readable by matlab
            disp(sprintf('couldnt read file %s as dicom', Filename_all{indextable(i,j)}));
        end
    end
end
dcm_vol = squeeze(dcm_vol);
volsize = size(dcm_vol);
cd(thisdir);