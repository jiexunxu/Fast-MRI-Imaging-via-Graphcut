function write_img(data, flstr, prec)
% by A Raj:  writes image block to binary file
% usage: write_img(data, flstr, prec); prec can be 'uint8', 'float32' etc

fid = fopen(flstr, 'wb');
fwrite(fid, data, prec);
fclose(fid);