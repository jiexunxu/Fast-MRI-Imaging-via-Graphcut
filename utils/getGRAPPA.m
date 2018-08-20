function x=getGRAPPA(GRAPPAPath)
    listing=dir(GRAPPAPath);
    pattern=listing(3).name;
    pattern=pattern(1:13);
    fname=strcat(GRAPPAPath, pattern, num2str(1), '.IMA');
    slice=dicomread(fname);
    x=zeros(size(slice, 1), size(slice, 2), length(listing));
    for i=1:length(listing)-2
        fname=strcat(GRAPPAPath, pattern, num2str(i), '.IMA');
        x(:, :, i)=dicomread(fname);
    end
    x=x(:, size(x, 2):-1:1, size(x, 3):-1:1);
    x=permute(x, [2 3 1]);
end