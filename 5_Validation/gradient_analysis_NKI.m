%% calculate gradients

%construct individual FC and gradient
for i=1:length(list)
name1=list(i,1).name;
funcfile=['...\NKI\FunImgARWSDCF_resliced\' name1, '\rFiltered_4DVolume.nii'];
M=x_gen_matrix_voxel(maskfile,funcfile);
zFC=fisherR2Z(M);
N = connectivity2normangle(zFC);
	[emb,res] = mica_diffusionEmbedding(N);
    gradient.emb=emb;
    gradient.res=res;
	gradient_filename = ['...\NKI\individual_gradient','\','g',list(s).name];
	save(gradient_filename,'gradient');
	clear M zFC N emb res;
end

%individual gradient pipeline
install_neuro_tools
sublist=dir('sub*');
num_part = length(sublist);
Opt.mode = 'qsub';
Opt.max_queued = 50;
Opt.flag_pause = false;
Opt.path_logs='.../results/logs/';
for n=1:num_part
        Job_name = ['gradient',num2str(n)];
    pipeline.(Job_name).command = 'individual_gradient(opt.sublist,opt.n)';
    pipeline.(Job_name).opt.sublist = sublist;
    pipeline.(Job_name).opt.n = n;

end
psom_run_pipeline(pipeline,Opt);

% Align individual's gradient components 
cd('...\NKI\results\indi_gradient\');
sublist = dir('gc*');
for s = 1 : length(sublist)
    load(sublist(s).name);
    embeddings{s,:} = gradient.emb; % 
    clear gradient;
end
nIterations = 100;
[realigned, target, xfms] = mica_iterativeAlignment(embeddings,nIterations);
realigned_gradient.realigned = realigned; % 
realigned_gradient.xfms = xfms;
save realigned_gradient realigned_gradient;

%% corr gradient AND age
load('...\NKI\code\table_model_resort.mat');
load('...\NKI\results\realigned_gradient.mat')

%calculate range,mean,std,explanation ratio and dispersion
indi_grad1(:,:)=realigned_gradient.realigned(:,1,:)';
indi_grad2(:,:)=realigned_gradient.realigned(:,2,:)';
save('indi_grad1','indi_grad1');
save('indi_grad2','indi_grad2');

range_indi_grad(:,1)=range(indi_grad1,2);
range_indi_grad(:,2)=range(indi_grad2,2);
save('range_indi_grad','range_indi_grad');
mean_indi_grad(:,1)=mean(indi_grad1,1);
mean_indi_grad(:,2)=mean(indi_grad2,1);
save('mean_indi_grad','mean_indi_grad');
std_indi_grad(:,1)=std(indi_grad1,0,2);
std_indi_grad(:,2)=std(indi_grad2,0,2);
save('std_indi_grad','std_indi_grad');

gradient_order_resort=zeros(390,11);
for i=1:390
    A=realigned_gradient.xfms{i};
    gradient_order=abs(A);
    [max_gradient_order,index]= max(gradient_order);
    gradient_order_resort(i,:)=index; 
end
gradient_indi_explan=zeros(390,136);
for i=1:390
load(['...\NKI\results\indi_gradient\g',sublist{i}]);
explan=gradient.res.lambdas./sum(gradient.res.lambdas);
gradient_indi_explan(i,:)=explan;
end
gradient_indi_explan_resort=zeros(390,2); 
for i=1:390
gradient_indi_explan_resort(i,1)=gradient_indi_explan(i,gradient_order_resort(i,1));
gradient_indi_explan_resort(i,2)=gradient_indi_explan(i,gradient_order_resort(i,2));
end
save('gradient_indi_explan','gradient_indi_explan');
save('gradient_order_resort','gradient_order_resort');
save('gradient_indi_explan_resort','gradient_indi_explan_resort');

%mixed linear model
load('...\NKI\code\table_model_resort.mat');
parfor i=1:17673 
prediction=indi_grad2(:,i);
[age_t,age_pValue,age_beta,age_lme,age_dfe] = mixed_model(prediction,table_model_resort);
stresults{i}.b=age_beta;  
stresults{i}.pval=age_pValue;
stresults{i}.stats.dfe=age_dfe;
stresults{i}.stats.yr=residuals(age_lme);
stresults{i}.tstat=age_t;
t_test(i,1)=age_t;
t_test(i,2)=age_pValue;
end
save('age_related_indi_grad2_stresults','stresults');
save('t_test','t_test');

%visualize
mask_hdr = spm_vol(maskfile);
mask_vol = spm_read_vols(mask_hdr);
mask_ind = reshape(mask_vol>0,1,[]);
grad=zeros(1,length(mask_ind)); 
grad(mask_ind)=t_test(:,1);  
[dim1,dim2,dim3]=size(mask_vol);
grad_nii=reshape(grad,dim1,dim2,dim3);
mask_hdr.fname='...\age_indi_grad2_t.nii';
mask_hdr.dt=[16,0]; 
spm_write_vol(mask_hdr,grad_nii);


