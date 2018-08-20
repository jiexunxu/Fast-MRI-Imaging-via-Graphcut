function newmat = img_auto_translate_mat(coT1, voxel_sizes, outputfileprefix)

%The array coT1 points to the images that require shifting dcm2nii converted images to match T1
%template. RunAtSize gives option to normalize and atlas the image by specifying the size of the atlas

% @ E. LoCastro, J. Xu

tsp=[-0.5; -0.6872; -0.4066];
dimz=size(coT1);

if size(dimz) < 4
    newmat=zeros([4 4]);
else
    newmat=zeros([4 4 dimz(4)]);
end


for t=1:size(newmat,3)
%     for i=1:3
%         newmat(i,4,t)=tsp(i)*dimz(i)*voxel_sizes(i);
%         newmat(i,i,t)=voxel_sizes(i);
%     end
    newmat(1:3, 4, t)=tsp'.*dimz.*voxel_sizes;
    newmat(1:3, 1:3, t)=diag(voxel_sizes);
end

newmat(4,4,:)=1;

if nargin > 2
    for t=1:size(newmat,3)
        if t < 10
            if t == 1
                ordnum = '';
            else
                ordnum = ['0' num2str(t)];
            end
        else
            ordnum=num2str(t);
        end
        outputfilename=[outputfileprefix ordnum '.img'];
        V=struct('fname', outputfilename, 'dim', size(coT1(:,:,:,t)), 'mat', newmat(:,:,t), 'pinfo', [1 0 0]', 'dt', [spm_type('float32') 0], 'n', [1 1], 'descrip', '');
        spm_write_vol(V,coT1(:,:,:,t));
    end
end        

end
