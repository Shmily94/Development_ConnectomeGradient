% Dependency:
% 1. "micaopen-master/diffusion_map_embedding/"

%% Prepar calculate voxel_wise FC matrix
maskfile=('...\grey_mask_0.2_4mm_AAL90_mask.nii');
subname=dir('...\FunImgARWSDCF_Resliced');
subname(1:2)=[];
size=length(subname);
for i=2:length(subname)
name1=subname(i,1).name;
funcfile=['...\FunImgARWSDCF_Resliced\' name1, '\rFiltered_4DVolume.nii'];
M=x_gen_matrix_voxel(maskfile,funcfile);
zFC=fisherR2Z(M);
save(name1,'zFC','-v7.3');
end

%% Estimate the individual-level gradient components
cd(subpath);
for s = 1 : length(sublist)
    data=cell2mat(struct2cell(load(sublist(s).name)));
	N = connectivity2normangle(data);
	[emb,res] = mica_diffusionEmbedding(N);
    gradient.emb=emb;
    gradient.res=res;
	gradient_filename = ['...\indi_gradient','\','g',sublist(s).name];
	save(gradient_filename,'gradient');
	clear data N emb res;
end

%% Aligment individual's gradient components
% put individual's embeddings in cells
cd('...\indi_gradient\');
sublist = dir('g*');
for s = 1 : length(sublist)
    load(sublist(s).name);
    embeddings{s,:} = gradient.emb; 
    clear gradient;
end
nIterations = 100;
[realigned, xfms] = mica_iterativeAlignment(embeddings,nIterations);
realigned_gradient.realigned = realigned; 
realigned_gradient.xfms = xfms;
save realigned_gradient realigned_gradient;

%% visualization
maskfile=('...\grey_mask_0.2_4mm_AAL90_mask.nii');
mask_hdr = spm_vol(maskfile);
mask_vol = spm_read_vols(mask_hdr);
mask_ind = reshape(mask_vol>0,1,[]);
grad=zeros(1,length(mask_ind)); 
grad(mask_ind)=gradient.emb(:,1);
[dim1,dim2,dim3]=size(mask_vol);
grad_nii=reshape(grad,dim1,dim2,dim3);
mask_hdr.fname='...\xxx.nii';
mask_hdr.dt=[16,0]; 
spm_write_vol(mask_hdr,grad_nii);

