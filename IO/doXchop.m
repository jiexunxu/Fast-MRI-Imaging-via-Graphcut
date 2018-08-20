function Chopped = doXChop(Mat)
% function Chopped = doXChop(Mat)
% multiply every second column by -1
% I don't know why this is called X chop but it is consistent
% with Richard's usage
[rows,cols] = size(Mat);
Chopped = Mat;
i = 2:2:cols;
Chopped(:,i) = Chopped(:,i)*-1;
