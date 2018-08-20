function [NiftiFileName, NiftiPath] = save_dicom_series_to_analyze(dicom_filenames)
% converts dicom series to analyze format (1 .img file and 1 .hdr file)
% saved in same folder as dicoms
% check to make sure spm2 or spm5 are in matlab path
thisdir = pwd;
[dirname, nm, ext] = fileparts(dicom_filenames{1});
cd(dirname);
% img_before = get_niftifiles(dirname);
hdr = spm_dicom_headers(char(dicom_filenames));
spm_dicom_convert(hdr, 'all', 'flat', 'nii');
% img_after = get_niftifiles(dirname);
% % look for changes in files list
% NiftiFileName = setdiff(img_after, img_before);

mkdir('tmp');
movefile('*.nii', 'tmp/');
cd('tmp');
NiftiFileName = get_niftifiles(pwd);
if isempty(NiftiFileName)
    disp('New nifti file(s) were not detected - please check! Either a nifti file of this name already exists or a new nifti file could not be saved');
elseif length(NiftiFileName)==1
    NiftiFileName = char(NiftiFileName);
end

if nargout>1
    NiftiPath = pwd;
end

cd(thisdir);


function niftifls = get_niftifiles(direc)
D = dir(direc);
imgind = find(~[D.isdir]);
D = D(imgind);
fls = {D.name};

t = 0; niftifls = [];
for i = 1:length(fls)
    if length(fls{i}) > 4
    if strcmp(fls{i}(end-3:end), '.img') || strcmp(fls{i}(end-3:end), '.nii')
        t = t+1;
        niftifls{t} = fls{i};
    end
    end
end
