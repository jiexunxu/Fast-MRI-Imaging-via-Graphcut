function q = read_img(flstr, SIZ, prec)
% By A Raj: reads image block from binary file
% usage: data = read_img(flstr, siz, prec, dispflg)
% prec is 'uint8', 'float32' etc; 
% can read 2d, 3d or 4d images, depending on length(siz)
% if length(siz)==1, then reads [siz x siz] image

try
    q = imread(flstr); 
catch
    def_ans = cell(1,5);
    prompt = cell(1,5);
    if nargin==0, flstr = '';  end
    title = 'Enter image attributes';
    prompt = {'File Name', 'rows', 'cols', 'slices', 'phases', 'Precision'};
    def_ans = {flstr, '0', '0', '0', '0', 'int32'};
    if nargin<2        
        ans = inputdlg(prompt, title, 1, def_ans, 'on');
        flstr = ans{1};
        nrows = str2num(ans{2});
        ncols = str2num(ans{3});
        nslices = str2num(ans{4});
        nphases = str2num(ans{5});
        prec = ans{5};
        if nphases==0
            if nslices==0
                SIZ = [nrows, ncols];
            else
                SIZ = [nrows, ncols, nslices];
            end
        else
            SIZ = [nrows, ncols, max(nslices,1), nphases];            
        end
    elseif nargin<3
       prec = char(inputdlg({'Enter image precision'}, title, 1, {'uint8'}));
    end
    fid = fopen(flstr);
    q = fread(fid, SIZ, prec);
    fclose(fid);    
end



