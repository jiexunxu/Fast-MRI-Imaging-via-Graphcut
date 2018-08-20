function [spmLabelMap, x, segMap]=spmVolumetrics2(x)
% SPM assumes that the brain data runs from bottom of head to top of head in z
% direction, back of head to front of head in y direction, and left to right of head in
% x direction.

% SIMMENS data runs from top to bottom of head in z direction, back of head to
% front in x direction, and left to right in y direction

% To correct orientation, switch x and y axis, and reverse z axis
%    x=permute(x, [2 1 3]);x=x(:, :, size(x, 3):-1:1);
%    x=double(abs(x));
    img_auto_translate_mat(x,[1 1 1], [pwd filesep 'test']);
    [spmLabelMap, segMap]=RunAtlasing_SVN([pwd filesep 'test.img'], 116);
%    rmdir('Segmented', 's');
   % rmdir('Normalized', 's');
 %   rmdir('Atlased116', 's');
end

function img_auto_translate_mat(coT1, voxel_sizes, outputfileprefix)
    tsp=[-0.5; -0.6872; -0.4066];
    dimz=size(coT1); 
    newmat=zeros([4 4]);
    newmat(1:3, 4)=tsp'.*dimz.*voxel_sizes;
    newmat(1:3, 1:3)=diag(voxel_sizes);
    newmat(4,4,:)=1;
    outputfilename=[outputfileprefix '' '.img'];
    V=struct('fname', outputfilename, 'dim', size(coT1), 'mat', newmat, 'pinfo', [1 0 0]', 'dt', [spm_type('float32') 0], 'n', [1 1], 'descrip', '');
    spm_write_vol(V,coT1);
end

function [spmLabelMap, segMap]=RunAtlasing_SVN(T1FileName, atlassize)
    %This is the directory of the toolbox.    
    SPMdir = 'X:\Projects\ApproximateNullVec\deployment\libs\spm8_pipeline'; 
    
    %The atlas to be used in the parcellation of the GM.    
    AtlasFileName = 'X:\Projects\ApproximateNullVec\deployment\libs\spm8_pipeline\atlas116.nii';

    if strcmpi(computer, 'GLNXA64')
        SPMdir = '/home/jiexunxu/SubspaceGraphcut/deployment/libs/spm8VolumetricsFiles'; 
        AtlasFileName = '/home/jiexunxu/SubspaceGraphcut/deployment/libs/spm8VolumetricsFiles/atlas116.nii';
    end
    %Leave this variable empty so that the subdirectories and output files are 
    %creaed in the patient's T1 directory. 
    %StrOutputDir = ['Atlasing' num2str(atlassize)];
    StrOutputDir = '';
    %This is the treshold for the GM segmentation (usually leave empty).
    Thresh = [];
    %Input the T1FileName (with directory) and produce the normalized T1, as
    %well as the WM, GM, CSF segmentation and parcellated GM atlas. Note: 
    %possible to pass cell array of T1filedies at once
    [AtlasedFileName, segMap] = Atlasing_spm8_SVN(T1FileName, StrOutputDir, AtlasFileName, 1, Thresh,SPMdir,0,atlassize,1);
    %Take each of the atlased T1's and count the number of voxels that have
    %been assigned to each of the 116 regions.
    spmLabelMap = spm_read_vols(spm_vol(deblank(AtlasedFileName)));   
    segMap=segMap+spmLabelMap;
end

function [GMSegFile, WMSegFile, CSSegFile, AtlasedFile] = Atlasing_spm8_SVN(Images, Output_dir, AtlasFile, TempFile, Thresh, SPMdir, MRIOPT,atlassize,affinereg)
%
% Syntax :
% atlasing(Images, Output_dir, AtlasFile, TempFile, Thresh);
%
% This script was developed over SPM2 toolbox and it computes automatically
% individual atlases based on Magnetic Resonance Images.
% The first step is based on the MRI image normalization to a stereotaxic
% space, MNI space(Montreal Neurological Institute). Here a transformations
% matrix is obtained. Then the MRI individual files are segmented in three
% different brain tissues (cerebral spinal fluid, gray and white matter) at
% this stage. During the second step each gray matter voxel is labeled with
% one structure label using the transformation matrix obtained in the
% normalization process and an anatomical atlas constructed by manual segmentation
% for a group of subjects.
%
% Input Parameters:
%   Images     : Individual MRI files.
%   Output_dir : Output directory for segmented, normalized and atlased
%                files. If the user doesn't change the output directory,
%                the resulting files are saved in the same address than the
%                individual MRI files.
%   AtlasFile  : Reference Atlas File used in the automatic labelling step.
%   TempFile   : Template file for the normalisation step. The user can use
%                any of the templates included at the SPM package or select
%                another, always taking care that this template file is in
%                the same space than the atlas that will be used in the
%                labelling step.
%   Thresh     : Threshold for gray matter segmentation.
%                Just the voxels with  higher probability than the threshold
%                are taken into acount in the automatic labelling step. If
%                the threshold isn't specified then an automatic one is taken.
%                All voxels with higher gray matter probabillity than 1-(GM+WM+CSF)
%                are taken into account in the automatic labelling step.
%                Being:
%                GM(A voxel V belongs to Gray Matter tissue with a probability GM)
%                WM(A voxel V belongs to White Matter tissue with a probability WM)
%                CSF(A voxel V belongs to Cerebral Spinal Fluid with a probability CSF)
%
% Related references:
% 1.- Ashburner J, Friston K. Multimodal image coregistration and partitioning--
%     a unified framework. Neuroimage. 1997 Oct;6(3):209-17.
% 2.- Voxel-based morphometry--the methods. Neuroimage. 2000 Jun;11(6 Pt
% 1):805-21.
% 3.- Evans AC, Collins DL, Milner B (1992). An MRI-based Stereotactic Brain
%     Atlas from 300 Young Normal Subjects, in: Proceedings of the 22nd Symposium
%     of the Society for Neuroscience, Anaheim, 408.
%
% See also: spm_normalise  spm_segment  Auto_Labelling
%__________________________________________________
% Authors: Lester Melie Garcia & Yasser Alem?? G??ez
% Neuroimaging Department
% Cuban Neuroscience Center
% November 15th 2005
% Version $1.0

%startup_varsonly;
if nargin<7, MRIOPT=0; end

spm_get_defaults;
warning off;
global defaults;
dseg = defaults.preproc;
dseg.output.GM = [0 0 1];
dseg.output.WM = [0 0 1];
dseg.output.CSF = [0 0 1];
dseg.tpm   = char(...
    fullfile(SPMdir,'tpm','grey.nii'),...
    fullfile(SPMdir,'tpm','white.nii'),...
    fullfile(SPMdir,'tpm','csf.nii'));
if affinereg == 0
    dseg.regtype = 'none';
end
%=====================Checking Input Parameters===========================%
if nargin==0
    [Images,~] = spm_select([1 Inf],'image','Selecting UnNormalized Images','',cd);
    [AtlasFileName,AtlasFilePath] = uigetfile({'*.img'},'Reading Reference Atlas File ...');
    AtlasFile = [AtlasFilePath AtlasFileName];
    [TempFile] = Choosing_TempFile;
    Thresh = input('Please select a threshold for gray matter segmentation:   ');
else
    if isempty(Images)
        [Images,~] = spm_select([1 Inf],'image','Selecting UnNormalized Images','',cd);
    end
    if isempty(AtlasFile)
        [AtlasFileName,AtlasFilePath] = uigetfile('*.img','Reading Reference Atlas File ...');
        AtlasFile = [AtlasFilePath AtlasFileName];
    end
    if isempty(TempFile)
        [TempFile] = Choosing_TempFile;
    end
    if ~exist('Thresh', 'var')
        Thresh = input('Please select a threshold for gray matter segmentation:   ');
    end
end
%=========================================================================%
%
%=========================Main program=====================================

V = spm_vol(Images);
Output_dir = char(Output_dir);
Ns = length(V);

%Preallocate stuff
GMSegFile=cellstr(char(zeros(Ns,1)));
WMSegFile=GMSegFile;
CSSegFile=GMSegFile;
AtlasedFile=GMSegFile;

for i=1:Ns
    disp(['Case ---> ' num2str(i)]);
    [pth, fn, xt] = fileparts(V(i).fname);
    ext=xt;
    subdir=pth;
    %Set the Atlasing output directory for Auto_Labelling function
    %The incoming variable 'Output_dir' now refers to a parent directory
    %housing all 3 (Normalized, Segmented & Atlasing) folders generated by
    %this function - EL
    if (nargin<2)||(isempty(Output_dir))
        Atlas_Output_dir = [pth filesep 'Atlased' num2str(atlassize)];
        mkdir(Atlas_Output_dir);
    else
        normfileout=[subdir filesep Output_dir];
        Atlas_Output_dir=[normfileout filesep 'Atlased' num2str(atlassize)];
        mkdir(normfileout);
        mkdir(Atlas_Output_dir);
        pth=normfileout;
    end
    %%%%%%%%%--------------- Normalization & Segmentation----------------%%%%%%%%%%%%%%%%
   
    disp('Normalizing and Segmenting ...');    
    mkdir(pth,'Normalized');
    mkdir(pth,'Segmented');
    V0 = spm_preproc(V(i),dseg);
    
    % convert segmentation information to sn-files
    [po,pin] = spm_prep2sn(V0);
        
    opts2.biascor   = 0;         %write bias corrected images    
    opts2.GM        = [0 0 1];   %write modulated, unmodulated normalized, and native GM
    opts2.WM        = [0 0 1];   %write modulated, unmodulated normalized, and native WM
    opts2.CSF       = [0 0 1];   %write modulated, unmodulated normalized, and native CSF
    opts2.cleanup   = 0;         %Do not cleanup

    % write segmented data
    spm_preproc_write(po,opts2);
    
    %save the segmentation information
    
    save(fullfile([pth filesep 'Normalized'], [fn '_seg_sn.mat']), '-struct', 'po');
    save(fullfile([pth filesep 'Normalized'], [fn '_seg_inv_sn.mat']), '-struct', 'pin');
    
    %Write the normalized T1 image to the subdirectory. 
    spm_write_sn(V(i).fname,po);
    if strcmp(xt, '.img')
        movefile([subdir filesep 'w' fn ext],[pth filesep 'Normalized']);
        movefile([subdir filesep 'w' fn '.hdr'],[pth filesep 'Normalized']);
    elseif strcmp(xt, '.nii')
        movefile([subdir filesep 'w' fn ext],[pth filesep 'Normalized']);
    end
    
    
    %Save the segmented images to the subdirectory.
    if strcmp(ext,'.img')
        movefile([subdir filesep 'c1' fn ext],[pth filesep 'Segmented']);
        movefile([subdir filesep 'c1' fn '.hdr'],[pth filesep 'Segmented']);
        movefile([subdir filesep 'c2' fn ext],[pth filesep 'Segmented']);
        movefile([subdir filesep 'c2' fn '.hdr'],[pth filesep 'Segmented']);
        movefile([subdir filesep 'c3' fn ext],[pth filesep 'Segmented']);
        movefile([subdir filesep 'c3' fn '.hdr'],[pth filesep 'Segmented']);
    elseif strcmp(xt, '.nii')
        movefile([subdir filesep 'c1' fn ext],[pth filesep 'Segmented']);
        movefile([subdir filesep 'c2' fn ext],[pth filesep 'Segmented']);
        movefile([subdir filesep 'c3' fn ext],[pth filesep 'Segmented']);
    end
    %define the images to normalize
    grey = fullfile([pth filesep 'Segmented'], ['c1' fn ext]);
    white = fullfile([pth filesep 'Segmented'], ['c2' fn ext]);
    csf = fullfile([pth filesep 'Segmented'], ['c3' fn ext]);
    Vgrey = spm_vol(grey);
    Vwhite = spm_vol(white);
    Vcsf = spm_vol(csf);
    
    %write unmodulated normalized
    spm_write_sn(Vgrey,po);
    spm_write_sn(Vwhite,po);
    spm_write_sn(Vcsf,po);
    
        %Open normalized image in mricron for quality check
    if MRIOPT
        openMRIcron([pth filesep 'Normalized' filesep 'w' fn ext ' -b 20 -o ' pth filesep 'Segmented' filesep 'wc1' fn ext]);
    end

    %%%%%%%%%--------- End of Normalization/Segmentation -----------%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%---------------  Atlasing --------------------%%%%%%%%%%%%%%%%
    disp('Atlasing ...');
    GMSegFile{i} = [pth filesep 'Segmented' filesep 'c1' fn ext];
    WMSegFile{i} = [pth filesep 'Segmented' filesep 'c2' fn ext];
    CSSegFile{i} = [pth filesep 'Segmented' filesep 'c3' fn ext];
    Transf_matname = fullfile([pth filesep 'Normalized'], [fn '_seg_sn.mat']);
    AtlasedFile{i} = Auto_Labelling(GMSegFile{i}, WMSegFile{i}, CSSegFile{i}, AtlasFile, Transf_matname, Atlas_Output_dir,Thresh);
    %%%%%%%%%--------- End of the Atlasing Step -----------%%%%%%%%%%%%%%%%
end;

GMSegFile = char(GMSegFile);
WMSegFile = char(WMSegFile);
CSSegFile = char(CSSegFile);
AtlasedFile = char(AtlasedFile);
%========================End of main program==============================%
return;

% replaced internal fn below choosing... by m-file... - AR
end

function [AtlasedFileName, VAtlas] = Auto_Labelling(GMSegFile, WMSegFile, CSSegFile, AtlasFile, Transf_matname, Output_dir, Thresh)
%
% Syntax : 
% Atlas =   Auto_Labelling(GMSegFile, WMSegFile, CSSegFile, AtlasFile,
% transf_matname)
%
% This function returns an individual Atlas in Analyze format. 
% The first step is to find the voxels that belongs to gray matter tissue 
% using the tissues probabilities maps previously obatained by the segmentation
% process. Due to the thresholding process, some holes, as well as isolated 
% points are present at gray matter volume. To solve this problem is used the 
% matlab function "imfill" to refills the internal holes and an internal
% function to reduce the isolated points(Uncoment if you want to employ it).
% Then each gray matter voxel is labeled with one structure label using the
% space transformation matrix (obtained in the normalization step), and an 
% anatomical atlas (constructed by manual segmentation for a group of
% subjects).
%
% Input Parameters:
%  GMSegFile      : Gray Matter Segmentation File.
%  WMSegFile      : White Matter Segmentation File.
%  CSSegFile      : Cerebral Spinal Fluid Segmentation File.
%  AtlasFile      : Reference Atlas File.
%  Transf_matname : Normalization Transform File.
%  Output_dir     : Output directory for segmented, normalized and atlased
%                    files. If the user doesn't change the output directory,
%                    the resulting files are saved in the same address than the
%                    Gray Matter Segmentation File.
%  Thresh         : Threshold for gray matter segmentation. 
%                   Just the voxels with  higher probability than the threshold
%                   are taken into acount in the automatic labelling step. If 
%                   the threshold isn't specified then an automatic one is taken. 
%                   All voxels with higher gray matter probabillity than 1-(GM+WM+CSF)
%                   are taken into account in the automatic labelling step.
%                   Being:
%                   GM(A voxel V belongs to gray matter tissue with a probability GM).
%                   WM(A voxel V belongs to white matter tissue with a probability WM). 
%                   CSF(A voxel V belongs to cerebral spinal fluid with a probability CSF).
% Output Parameter:
%   Atlas : Individual Gray Matter Atlas.
%
% Related References:
% 1.- Ashburner J, Friston K. Multimodal image coregistration and partitioning-- 
%     a unified framework. Neuroimage. 1997 Oct;6(3):209-17. 
% 3.- Evans AC, Collins DL, Milner B (1992). An MRI-based Stereotactic Brain
%     Atlas from 300 Young Normal Subjects, in: Proceedings of the 22nd Symposium 
%     of the Society for Neuroscience, Anaheim, 408.
%
% Note: The morphological treatment(erosion and dilation) are comented due 
%       to its dependence of the results with the structure element used to 
%       erode and dilate.      
%
% See also: imfill  imerode  imdilate  spm_segment atlasing spm_normalise
%__________________________________________________________________________
% Authors:  Yasser Alem?? G??ez & Lester Melie Garc??  
% Neuroimaging Department
% Cuban Neuroscience Center
% Last update: November 15th 2005
% Version $1.0
%=========================Main program=====================================
spm_defaults;
global defaults
defaults.analyze.flip = 0; % Not flipped Images

VSeg =  spm_vol(GMSegFile);
VWM = spm_vol(WMSegFile);
VCS = spm_vol(CSSegFile);
[Apth,Aname,Aext] = fileparts(AtlasFile);  % Here the matrix Iatlas is loaded.
Vatlas = spm_vol(AtlasFile);Iatlas = spm_read_vols(Vatlas); Atlastype = Aext(2:end);

%%%%%%%%%%%-----------------Thresholding-------------------%%%%%%%%%%%%%%%%
 Iseg_mask = zeros(VSeg.dim(1),VSeg.dim(2),VSeg.dim(3));
 AllBrainMask = Iseg_mask;
 eThresh = exist('Thresh', 'var');
 for z = 1:VSeg.dim(3)
     %% tissues
     GM  = spm_slice_vol(VSeg,spm_matrix([0 0 z]),VSeg.dim(1:2),0);
     WM  = spm_slice_vol(VWM,spm_matrix([0 0 z]),VWM.dim(1:2),0);
     CSF = spm_slice_vol(VCS,spm_matrix([0 0 z]),VCS.dim(1:2),0);
     None = 1 - (GM + WM + CSF);
     %% inds
     %Finds the voxel indices with highest probability of being Gray Matter
%      Thresh2 = Thresh;
%      size(GM); size(WM); size(CSF); size(None); size(Thresh);
     if (eThresh==0) || (isempty(Thresh)) || (Thresh == 0)
         indGM = find((GM > WM) & (GM > CSF) & (GM > None));
     else
         indGM  = find((GM>CSF)&(GM>WM)&(GM>=Thresh));
     end
     %%%% --------------- Creating WM? Why? Brain Mask  ------------------%%%%
     %%%%  Lester Melie  April 11th 2006
     indWM = find((WM > GM) & (WM > CSF) & (WM > None));  % Indices with highest probability of being White Matter
     %indCSF = find((CSF > GM) & (CSF > WM) & (CSF > None));
     if (~isempty(indGM))&&(~isempty(indWM)) %&(~isempty(indCSF))
         %WM(unique([indGM; indWM; indCSF])) = 1;
         WM(unique([indGM; indWM])) = 1; %Masks out both GM and WM voxel indices
         AllBrainMask(:,:,z) = WM; %Doesn't look like this is used again in this routine
     end;
     %% brain
     Brain = sparse(VSeg.dim(1),VSeg.dim(2));
     Brain(indGM) = 1; %GM mask from GM voxel indicesfrom:amy.kuceyeski@gmail.com
     Iseg_mask(:,:,z) = full(Brain);
end
%%%%%%%%%-------------- End of Thresholding ---------------%%%%%%%%%%%%%%%%
disp(['Refilling...']);  
Iseg_mask = imfill(logical(Iseg_mask),'holes');
clear IWM ICS Iseg ind WM CSF None;
%Iseg_mask = Iso_Rem(Iseg_mask,7);
ref=tic;VAtlas = Get_Norm_Coord(Atlastype,AtlasFile,Transf_matname,Iseg_mask, Iatlas, Output_dir);toc(ref);
AtlasedFileName = VAtlas.fname;
%========================End of main program==============================%
end

%=====================Internal functions==================================%
function Vol = Get_Norm_Coord(Atlastype, AtlasFile, Matname, Iseg_mask, Iatlas, Output_dir)
%
%
% Input Parameters:
%   Atlastype   : Gray Matter Segmentation File
%   AtlasFile   : Reference Atlas File
%   matname     : Normalisation Transform File
%   Iseg_mask   : White Matter Segmentation File <-- GM only?
%   Iatlas      : Cerebral Spinal Fluid Segmentation File <-- Atlas File
%   Output_dir  : Output Directory
%
% Note: This is based on the function 'spm_write_sn' developed by PhD.John Ashburner(FIL,UCL).
%
%__________________________________________________________________________
% Authors:  Yasser Alem?? G??ez & Pedro Vald?? Hern??dez  
% Neuroimaging Department
% Cuban Neuroscience Center
% Last update: November 15th 2005
% Version $1.0
warning off
global defaults
defaults.analyze.flip = 0; % Not flipped Images
Vatlas = spm_vol(AtlasFile);
if length(Vatlas.dim)==4
    dt = [Vatlas.dim(4) 0];
elseif length(Vatlas.dim)==3
    dt = Vatlas.dt;
end
% VAtlas_mat = Vatlas.mat;
load('-mat',Matname); 
[Vpth,Vname,Vext] = fileparts(VF.fname); mat = VF.mat; 
if strcmp(spm('ver'),'SPM2')
dim = [VF.dim(1:3) dt(1)];
elseif strcmp(spm('ver'), 'SPM5') || strcmp(spm('ver'), 'SPM8')
    dim = [VF.dim(1:3)];
end
Vol = struct('fname','','dim',dim,'mat',mat,'pinfo',[1 0 0]',...
    'descrip','Atlas image','dt',dt);
Vol.fname =[char(Output_dir) filesep Vname '_Atlas' Vext];
Vol = spm_create_vol(Vol);
x = 1:VG(1).dim(1); y = 1:VG(1).dim(2); z = 1:VG(1).dim(3);
if ~isempty(Tr)
    BX = spm_dctmtx(VG(1).dim(1),size(Tr,1),x-1);
    BY = spm_dctmtx(VG(1).dim(2),size(Tr,2),y-1);
    BZ = spm_dctmtx(VG(1).dim(3),size(Tr,3),z-1);
end
[X,Y] = ndgrid(x,y); clear x y
y1 = single(0); y1(VG(1).dim(1),VG(1).dim(2),VG(1).dim(3)) = 0;
y2 = single(0); y2(VG(1).dim(1),VG(1).dim(2),VG(1).dim(3)) = 0;
y3 = single(0); y3(VG(1).dim(1),VG(1).dim(2),VG(1).dim(3)) = 0;
M = VG(1).mat;
for j=1:length(z);
    if ~isempty(Tr)
        X1 = X    + BX*get_2Dtrans(Tr(:,:,:,1),BZ,j)*BY';
        Y1 = Y    + BX*get_2Dtrans(Tr(:,:,:,2),BZ,j)*BY';
        Z1 = z(j) + BX*get_2Dtrans(Tr(:,:,:,3),BZ,j)*BY';
    else
        X1 = X; Y1 = Y; Z1 = z(j);
    end
    y1(:,:,j) = single(M(1,1)*X1 + M(1,2)*Y1 + M(1,3)*Z1 + M(1,4));
    y2(:,:,j) = single(M(2,1)*X1 + M(2,2)*Y1 + M(2,3)*Z1 + M(2,4));
    y3(:,:,j) = single(M(3,1)*X1 + M(3,2)*Y1 + M(3,3)*Z1 + M(3,4));
end
clear X1 Y1 Z1 X Y z
disp(['Inverting the deformations field...']); 
M = Affine/VG(1).mat; M(4,:) = [0 0 0 1];
[iy1,iy2,iy3] = spm_invdef(y1,y2,y3,VF.dim(1:3),M,VG(1).mat); 
clear y1 y2 y3
M = inv(Vatlas.mat);
for j = 1:VF.dim(3)
    A = zeros(VF.dim(1),VF.dim(2));
    %disp(['Slice ----> ' num2str(j)]);
    if sum(sum(Iseg_mask(:,:,j  )))~=0
        X2 = M(1,1)*double(iy1(:,:,j)) + M(1,2)*double(iy2(:,:,j)) + M(1,3)*double(iy3(:,:,j)) + M(1,4);
        Y2 = M(2,1)*double(iy1(:,:,j)) + M(2,2)*double(iy2(:,:,j)) + M(2,3)*double(iy3(:,:,j)) + M(2,4);
        Z2 = M(3,1)*double(iy1(:,:,j)) + M(3,2)*double(iy2(:,:,j)) + M(3,3)*double(iy3(:,:,j)) + M(3,4);
        A = spm_sample_vol(Iatlas,X2,Y2,Z2,0);
    end
    spm_write_plane(Vol,squeeze(A.*Iseg_mask(:,:,j)),j);
end
fclose all;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function T2 = get_2Dtrans(T3,B,j)
 d   = [size(T3) 1 1 1];
 tmp = reshape(T3,d(1)*d(2),d(3));
 T2  = reshape(tmp*B(j,:)',d(1),d(2));
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [arr_names,sz]=Dir2Arr(srcpth,extns,omit)
%Dir2Arr sends dir search filenames to matlab array variable (like its pre-
%decessor, DirToArr), but with extended functionality. srcpth, extns and 
%omit can be all be (cell) arrays of text. Returns text array with
%spaced-padding at end
%
% @ E. LoCastro
MSL=500;
if nargin > 2
    omit_arr=Dir2Arr(srcpth,omit);
    if isempty(omit_arr)
        omit_list={};
    else
        O=size(omit_arr,1);
        for o=1:O
            [~,f,e]=fileparts(omit_arr(o,:));
            omit_arr(o,:)=[f e blanks(MSL-length([f e]))];
        end
        omit_list=cellstr(omit_arr);
    end
else
    omit_list={};
end

arr_names=[];
if ~iscellstr(srcpth), srcpth=cellstr(srcpth); end
if ~iscellstr(extns), extns=cellstr(extns); end   
dupe_prevention={};
for p=1:length(srcpth)
    read_dir=srcpth{p};
    for e=1:length(extns)
        G=dir([read_dir filesep extns{e}]);
        for g=1:length(G)
            if ~ismember(G(g).name,omit_list) && ~strcmp(G(g).name(1),'.') && ~ismember(G(g).name,dupe_prevention)
                fname=[read_dir filesep G(g).name];
                dupe_prevention{end+1}=G(g).name;
                arr_names=[arr_names; fname blanks(MSL-length(fname))];
            end
        end    
    end
    %disp(p);
end 
sz=size(arr_names,1); %sz=S(1);
disp(['Dir2Arr: ' int2str(sz) ' entries.']);
end