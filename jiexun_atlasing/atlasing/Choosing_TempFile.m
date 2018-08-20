function [TempFile] = Choosing_TempFile(Type)
%=====================Used to be Internal function==================================%
if nargin==0
    Type = input('Choosing template(T1, T2, PET, SPECT, EPI, PD or Other):   ','s');
    Type = lower(Type);
end
switch Type
    case 't1'
        if exist(fullfile(spm('Dir'),'templates','T1.nii'))
            TempFile = fullfile(spm('Dir'),'templates','T1.nii');
        else
            TempFile = fullfile(spm('Dir') , 'templates' , 'T1.mnc');
        end
    case 't2'
        if exist(fullfile(spm('Dir'),'templates','T2.nii'))
            TempFile = fullfile(spm('Dir'),'templates','T2.nii');
        else
            TempFile = fullfile(spm('Dir') , 'templates' , 'T2.mnc');
        end
    case 'pet'
        if exist(fullfile(spm('Dir'),'templates','PET.nii'))
            TempFile = fullfile(spm('Dir'),'templates','PET.nii');
        else
            TempFile = fullfile(spm('Dir') , 'templates' , 'PET.mnc');
        end
    case 'spect'
        if exist(fullfile(spm('Dir'),'templates','SPECT.nii'))
            TempFile = fullfile(spm('Dir'),'templates','SPECT.nii');
        else
            TempFile = fullfile(spm('Dir') , 'templates' , 'SPECT.mnc');
        end
    case 'epi'
        if exist(fullfile(spm('Dir'),'templates','EPI.nii'))
            TempFile = fullfile(spm('Dir'),'templates','EPI.nii');
        else
            TempFile = fullfile(spm('Dir') , 'templates' , 'EPI.mnc');
        end
    case 'pd'
        if exist(fullfile(spm('Dir'),'templates','PD.nii'))
            TempFile = fullfile(spm('Dir'),'templates','PD.nii');
        else
            TempFile = fullfile(spm('Dir') , 'templates' , 'PD.mnc');
        end
    case 't1-canonical'
        TempFile='/media/94f7fc1e-fa27-4b24-99f3-4b461665a4a4/home/eve/Documents/MATLAB/eve_tools/mricron/templates/MNI152_T1_1mm.nii.gz';
    case 't1_brain'
        TempFile='/media/94f7fc1e-fa27-4b24-99f3-4b461665a4a4/home/eve/Documents/MATLAB/spm8/templates/MNI152lin_T1_1mm_brain.nii';
    case 'other'
        [TempFileName,TempFilePath] = uigetfile('*.mat;*.img','Reading Reference Template File ...');
        TempFile = [TempFilePath TempFileName];
    case 'infant'
        TempFile= fullfile(spm('Dir'),'templates','nihpd_sym_04.5-08.5_t1w.nii');
    otherwise
        disp('Unknown template file ...');
end;
return;

