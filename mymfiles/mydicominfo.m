function [info, file] = mydicominfo(msgname, varargin)
%DICOMINFO  Read metadata from DICOM message.
%   INFO = DICOMINFO(FILENAME) reads the metadata from the compliant
%   DICOM file specified in the string FILENAME.
%
%   INFO = DICOMINFO(FILENAME, 'dictionary', D) uses the data dictionary
%   file given in the string D to read the DICOM message.  The file in D
%   must be on the MATLAB search path.  The default value is dicom-dict.mat.
%
%   Example:
%
%     info = dicominfo('CT-MONO2-16-ankle.dcm');
%
%   See also DICOMREAD.

%   INFO = DICOMINFO(FILE) reads the metadata from the currently open
%   message in the FILE structure.  This syntax allows DICOMREAD to call
%   DICOMINFO without reopening the file.
%
%   [INFO, FILE] = DICOMINFO(FILE) also returns the file structure (FILE)
%   in case the warning information has changed.

%   Copyright 1993-2002 The MathWorks, Inc.
%   $Revision: 1.6 $  $Date: 2002/03/15 15:27:06 $


% This function (along with DICOMREAD) implements the M-READ service.
% This function implements the M-INQUIRE FILE service.


%
% Parse input arguments.
%

if (nargin < 1)
    
    error('DICOMINFO requires at least one argument.')
    
end

[args, msg] = parse_inputs(msgname, varargin{:});

if (~isempty(msg))
    error(msg);
end

%
% Determine what to do based on whether file struct was passed in.
%

if (isstruct(msgname))
    % Get info on an already opened message.
    
    file = msgname;
    
    % Reset number of warnings.
    file = dicom_warn('reset', file);

    %
    % Extract the metadata.
    %
    
    % Create an Info structure with default values.
    info = dicom_create_meta_struct(file);
    
    % Look for file metadata.
    if (dicom_has_fmeta(file))
        
        [info, file] = dicom_read_fmeta(file, info, args.Dictionary);
        
    end
    
    % Find the encoding for the remaining metadata.
    file = dicom_set_mmeta_encoding(file, info);
    
    % Extract the message metadata.
    [info, file] = dicom_read_mmeta(file, info, args.Dictionary);
    
    % Assign the IMFINFO specific values.
    info = dicom_set_imfinfo_values(info, file);
   
else
    % Message string was passed in by user.

    %
    % Create File structure with uninitialized values.
    %
    
    file = dicom_create_file_struct;
    
    file.Messages = msgname;
    file.Location = args.Location;

    % Reset number of warnings.
    file = dicom_warn('reset', file);

    %
    % Get message to read.
    %
    
    file = dicom_get_msg(file);
    
    if (isempty(file.Messages))
        
        if (isequal(file.Location, 'Local'))
            msg = sprintf('File "%s" not found.', msgname);
        else
            msg = 'Query returned no matches.';
        end
        
        error(msg)
        
    end
    
    % Create a container for the output.
    info = {};
    
    %
    % Read metadata from each message.
    %
    
    for p = 1:length(file.Messages)
        
        %
        % Open the message.
        %
        
        file = dicom_open_msg(file, 'r');
        
        %
        % Extract the metadata.
        %
        
        % Create an Info structure with default values.
        info{p} = dicom_create_meta_struct(file);
        
        % Look for file metadata.
        if (dicom_has_fmeta(file))
            
            [info{p}, file] = dicom_read_fmeta(file, info{p}, args.Dictionary);
            
        end
        
        % Find the encoding for the remaining metadata.
        file = dicom_set_mmeta_encoding(file, info{file.Current_Message});
        
        % Extract the message metadata.
        [info{p}, file] = dicom_read_mmeta(file, info{p}, args.Dictionary);
        
        % Assign the IMFINFO specific values.
        info{p} = dicom_set_imfinfo_values(info{p}, file);
   
        %
        % Close the message.
        %
        
        file = dicom_close_msg(file);
        
    end  % For
    
    % Remove from cell array if only one structure.
    if (file.Current_Message == 1)
        info = info{1};
    end
    
end  % (isstruct(msgname))



%%%
%%% Function parse_inputs
%%%
function [args, msg] = parse_inputs(msgname, varargin)

% Set default values
args.Dictionary = 'dicom-dict.mat';

msg = '';

% Determine if messages are local or network.
% Currently only local messages are supported.
args.Location = 'Local';

% Parse arguments based on their number.
if (nargin > 1)
    
    paramStrings = {'dictionary'};
    
    % For each pair
    for k = 1:2:length(varargin)
       param = lower(varargin{k});
       
            
       if (~ischar(param))
           msg = 'Parameter name must be a string';
           return
       end

       idx = strmatch(param, paramStrings);
       
       if (isempty(idx))
           msg = sprintf('Unrecognized parameter name "%s"', param);
           return
       elseif (length(idx) > 1)
           msg = sprintf('Ambiguous parameter name "%s"', param);
           return
       end
    
       switch (paramStrings{idx})
       case 'dictionary'

           if (k == length(varargin))
               msg = 'No data dictionary specified.';
               return
           else
               args.Dictionary = varargin{k + 1};
           end

       end  % switch
       
    end  % for
           
end

	

