function [I, J, K, vals] = get_user_ROI(dirname)
% Written by Ashihs Raj March 10 2010
% Returns ROI indices and intensity values
% Added DICOM read, manual slice select
% dirname is where the case resides
%   Detailed explanation goes here (maybe someday)

codedir = pwd;
cd(codedir);
if nargin < 1
    dirname = 'C:\coolcap\example_data';
end
thresh_ratio = 0.06;

% ask for input directory 
user_path = uigetdir(dirname, 'Select directory containing dicoms of interest (CANCEL to accept suggested)');
if ~isequal(user_path,0)    
    dirname = user_path;
end
cd(dirname);
[read_success, found_dicom] = preview_dicoms;
if ~read_success
    if found_dicom
        % ask again for input directory 
        user_path = uigetdir(dirname, 'Try again: Select directory containing dicoms of interest (CANCEL to accept suggested)');
        if ~isequal(user_path,0)    
            dirname = user_path;
        end
        cd(dirname);
    else
        error('selectedt folder does not contain subfolders with readable DICOMS!');
    end
end

cd(dirname);
D = dir(pwd);
imgind = find(~[D.isdir]);
info = dicominfo(D(imgind(1)).name);
if ~strcmp(info.SeriesDescription, 'AXL DWI')
    disp('WARNING: this does not appear to be Axial DWI data!!');
end
nrows = info.Rows;
ncols = info.Columns;
nslices = info.ImagesInAcquisition,
imagesize = [nrows, ncols, nslices];
outimages = zeros(imagesize);
outmasks = zeros(imagesize);
for i = 1:length(imgind)
    info = dicominfo(D(imgind(i)).name);
    imgnum(i) = info.InstanceNumber;
    rawdata(:,:,i) = single(dicomread(info));
end
rawdata = rawdata(:,:, imgnum);
[nrows, ncols, nslices] = size(rawdata);

mask=(rawdata>=thresh_ratio*max(abs(rawdata(:))));
rawdata=mask.*rawdata;
    
    % Now get user to input ROIs 
    figure; montage(reshape(uint8(255*rawdata/max(rawdata(:))), [imagesize(1:2), 1, imagesize(3)])); 
    title('images, left to right, top to bottom');
    h1 = gcf;
    vals = []; 
    inds = [];
    user_done = false; 
    while ~user_done
        user_slice = sscanf(char(inputdlg(sprintf('input slice number to select ROI'))), '%03d');
        q = rawdata(:,:,user_slice);
        figure; imshow(uint8(255*q/max(q(:)))); truesize(size(q)*1.5); h2 = gcf;   zoom(2); 
        figure(h2);
        [bw, xi, yi] = roipoly;
        inds = [inds; find(bw)];
        vals = [vals; q(find(bw))];
        button = questdlg('Select another ROI?','ROI Select', 'No');
        if strcmp(button, 'No') | strcmp(button, 'Cancel')
            user_done = true;
        end
    end
    delete(h1); delete(h2);
    [I, J, K] = ind2sub(size(rawdata), inds);