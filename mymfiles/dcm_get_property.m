function [val, err] = dcm_get_property(dcm_filename, names)
% DCM_GET_PROPERTY: Read a property out of a DICOM properties array.
% Raj added multimple fields search - pass name as string cell

if iscell(names), is_cell = 1; else is_cell = 0; end
[properties] = dcm_properties(dcm_filename);
names = cellstr(names);
err = 0;

for kk = 1:length(names)
    name = char(names{kk});
    found_prop = 0;
    lidx = 0;
    ridx = length(properties) + 1;

    % Perform a binary search to find the property's value.
    while abs(lidx - ridx) > 1
      idx = floor((ridx + lidx) / 2);
      order = sort({properties{idx}.name, name});

      properties{idx};
      if (strcmpi(properties{idx}.name, name))
        val{kk} = properties{idx}.value;
        found_prop = 1;
        break;
      elseif (strcmpi(order{1}, name) )
        ridx = idx;
      else
        lidx = idx;
      end
    end

    if found_prop == 0 % The property was not in the list.
        %disp(['Field ' name ' was not found using dcm_get_property']);
        val{kk} = '';
        err = 1;
    end
end

if is_cell == 0, val = char(val); end
