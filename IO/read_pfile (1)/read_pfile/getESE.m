function esever = getESE( pfile )
%GETESE Summary of this function goes here
%   Detailed explanation goes here

fid=fopen(pfile, 'rb');
rev=fread(fid,4,'*uchar');
fclose(fid);
rev = round(10*mytypecast(rev, 'float32'));
switch rev
    case 143
        esever='14M5';
    case 150
        esever='15M3';
    case 160
        esever='16M3';
    case 200
        esever='20M3';
    case 220
        esever='22M4';
    otherwise
        error('unknown ESE version');
end
end

