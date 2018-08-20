function Chopped = doYChop(Mat)
% function Chopped = doYChop(Mat)
% multiply every second row by -1
% I don't know why this is called Y chop but it is consistent
% with Richard's usage
[rows,cols] = size(Mat);
Chopped = Mat;
i = 2:2:rows;
Chopped(i,:) = Chopped(i,:)*-1;
   