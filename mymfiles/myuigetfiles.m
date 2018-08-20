function [FileNames, PathName] = myuigetfiles(filt_str, msg)

% tries to reproduce known functionality of uigetfiles - 
% cell if multiple, non-cell string if single
% also removes the wrapping effect 

% [user_datafile, user_datapath] = uigetfiles(filt_str, msg);
[user_datafile, user_datapath] = uigetfile(filt_str, msg, 'Multiselect', 'on');
if (isequal(user_datafile,0) || isequal(user_datapath,0))
    disp('no file(s) selected');
    FileNames = 0;
    PathName = 0;
else
    FileNames = user_datafile;
    PathName = user_datapath;    
    if iscell(user_datafile) 
        if length(user_datafile)>1
            listlen = length(user_datafile);
           % wrap list by one - needed for some reason...
           tmp = user_datafile;
           for kk = 1:listlen-1, FileNames{kk} = tmp{kk+1};  end
           FileNames{listlen} = tmp{1};
        else
            FileNames = char(user_datafile);
        end
    else
        FileNames = user_datafile;
    end
end

