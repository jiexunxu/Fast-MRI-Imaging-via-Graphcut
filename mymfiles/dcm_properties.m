function [properties] = dcm_properties(filename)
% DCM_PROPERTIES: Read fields out of a DICOM file's header.

% TODO: Set this variable to the path to your copy of dcmdump.
DCMDUMP='C:\dcmtk-3.5.4-win32-i386\bin\dcmdump.exe';

% Check to see if dcmdump is valid.
if (exist(DCMDUMP, 'file') ~= 2)
    error(['Unable to access dcmdump at ' DCMDUMP 10]);
end

% Call dcmdump.
[status, output] = system([DCMDUMP ' ' filename]);
if (status ~= 0)
    error(['dcmdump failed on file ' filename]);
end

% Use regular expressions to pull the field names and values out of the output
% from dcmdump.
lines = mysplit(output);
properties = {};
names = {};
for i = 1:length(lines)
  line = lines(i, :);
  vals = regexp(line, '\[(.*)\]\s+#\s+\d+, \d+ ([^ ]+)', 'tokens');
  if length(vals) > 0
    newIdx = length(properties) + 1;
    properties{newIdx}.name = vals{1}{2};
    properties{newIdx}.value = vals{1}{1};
    names{newIdx} = [vals{1}{2} sprintf('#%d', newIdx)];
  end
end

% Sort the properties by name (this is gross because Matlab's sort is stupid).
names = sort(names);
sortedProperties = {};
for i = 1:length(names)
  srcIdx = sscanf(names{i}(strfind(names{i}, '#') + 1:end), '%d');
  sortedProperties{i} = properties{srcIdx};
end
properties = sortedProperties;

