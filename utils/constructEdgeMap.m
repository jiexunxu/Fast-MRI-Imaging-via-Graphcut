function edgeMap=constructEdgeMap(x0, thres)
    [numRows, numCols, numSlices]=size(x0);    
    edgeMap=zeros(numRows, numCols, numSlices);
    x0=abs(x0);
    for i=1:numSlices
        bw=myCanny(x0(:, :, i), 'canny', thres, 0.05);  
      %  bw=myCanny(x0(:, :, i), 'canny', 0.0001, 0.05);        
      %  bw=edge(x0(:, :, i), 'canny', 0.01, 0.2);
        edgeMap(:, :, i)=bw;
    end    
  %  for i=1:numSlices
  %      edgeMap(:, :, i)=bwmorph(edgeMap(:, :, i), 'clean', 3);
  %      edgeMap(:, :, i)=bwmorph(edgeMap(:, :, i), 'spur', 3);
  %      edgeMap(:, :, i)=bwmorph(edgeMap(:, :, i), 'clean', 3);
  %  end   
end