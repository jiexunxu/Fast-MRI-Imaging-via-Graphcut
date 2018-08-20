function [dcm_vol, volsize, voxsize, ndicoms] = myread_dicoms(data_dirname)

thisdir = pwd;
% ask for input directory 
if nargin==0
    data_dirname = thisdir;
end
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
    [val, err] = dcm_get_property(D(imgind(i)).name, {'InstanceNumber', 'AcquisitionNumber', 'Filename', 'SliceThickness', 'PixelSpacing', 'ImagePositionPatient'});
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
                acq_all(i) = info.AcquisitionNumber;
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
nacq = max(acq_all);
if ~(length(imgind) == max(ins_all) || length(imgind) == max(ins_all) * nacq)
    disp('WARNING: discrepancy between noof files, InstanceNumber and AcquisitionNumber');
    nslices = floor(length(imgind)/nacq);
    nooffiles = length(imgind), mxInstanceNumber = max(ins_all), nacq, nslices,
else
    nslices = min(max(ins_all), floor(length(imgind)/nacq));    
end
q = single(dicomread(Filename_all{1}));
[nrows, ncols] = size(q);
dcm_vol = zeros(nrows, ncols, nslices, nacq, 'single');

for acq = 1:nacq
    ii = find(acq_all == acq);
    ins = ins_all(ii);
    [st, iii] = sort(ins);
    ii = ii(iii);
    flist = {Filename_all{ii}};
    for sl=1:length(ii)
        try
            q = single(dicomread(flist{sl}));
            dcm_vol(:,:,sl,acq) = q;
        catch
            % dont read this file - not dicom or not readable by matlab
            disp(sprintf('couldnt read file %s as dicom', flist{sl}));
        end
    end
end
dcm_vol = squeeze(dcm_vol);
volsize = size(dcm_vol);
cd(thisdir);