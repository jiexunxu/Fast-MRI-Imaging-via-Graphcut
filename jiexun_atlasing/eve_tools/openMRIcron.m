function openMRIcron(MRIimg, MNItemptype)
% Calls LINUX version of MRIcron. May need modification for Windows usage.
% http://www.mccauslandcenter.sc.edu/mricro/mricron/install.html
% @ E. LoCastro

startup_varsonly;

if nargin<2
    MNItemptype=0;
end

[d,~,e]=fileparts(MRIimg);

if isempty(d)
        MRIimg=[pwd filesep MRIimg];
end

if isempty(e)
        MRIimg=[MRIimg '.img'];
end

if MNItemptype==0
    for i=1:size(MRIimg,1)
        system([mricronpath filesep 'mricron ' deblank(MRIimg(i,:)) ' &']);
    end
else
    for i=1:size(MRIimg,1)
        system([mricronpath filesep 'mricron ' deblank(MRIimg(i,:)) ' -o ' Choosing_TempFile(MNItemptype) ' &']);
    end
end

end
