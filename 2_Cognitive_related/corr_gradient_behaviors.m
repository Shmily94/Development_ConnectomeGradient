%% corr working memory and primary-transmodal gradient measures

load('...\WM.mat');
load('...\gradient_WM.mat');

%exclude subjects and generate the table_model_resort_WM.
index1=WM(:,1)>0.4;
index2=WM(:,2)>0;
index3=WM(:,3)>0;
index=index1.*index2.*index3;
index=logical(index);
WM_exclude=WM(index,:);

gradient_WM_exclude=gradient_WM(index,:);
subname=table_model_resort_WM.subname(index);
sex=table_model_resort_WM.sex(index);
meanFD=table_model_resort_WM.meanFD(index);

table_model_resort_WM(length(WM_exclude)+1:end,:) = [];
table_model_resort_WM.subname=subname;
table_model_resort_WM.sex=sex;
table_model_resort_WM.meanFD=meanFD;
table_model_resort_WM.cov_age=[];
table_model_resort_WM.baseline=WM_exclude(:,7); %0-back dprime
save('table_model_resort_WM','table_model_resort_WM')
save('WM_exclude','WM_exclude');
save('gradient_WM_exclude','gradient_WM_exclude');

%corr WM and age
load('...\table_model_resort.mat');
table_model_resort.baseline=WM_exclude(:,7);
prediction=WM_exclude(:,9);
[age_t,age_pValue,age_beta,age_lme,age_dfe] = mixed_model_WM(age,prediction,table_model_resort);

%corr WM and global gradient topography
load('...\table_model_resort_WM.mat');
for j=1:9
for i=1:6
prediction=WM_exclude(:,j); 
grad=gradient_WM_exclude(:,i);
[grad_t,grad_pValue,grad_beta,grad_lme,grad_dfe] = mixed_model_WM(grad,prediction,table_model_resort_WM);
results_WM{j}(i,1)=grad_t;
results_WM{j}(i,2)=grad_pValue;
end
end

% corr WM and regional gradient scores
load('...\WM\index');
indi_grad2=(-1).* indi_grad2;
indi_grad2_exclude=indi_grad2(index,:);
prediction=WM_exclude(:,9); 
stresults=cell(17673,1);

parfor i=1:17673 
grad=indi_grad2_exclude(:,i);
[grad_t,grad_pValue,grad_beta,grad_lme,grad_dfe] = mixed_model_WM(grad,prediction,table_model_resort_WM);
stresults{i}.b=grad_beta;  
stresults{i}.pval=grad_pValue;
stresults{i}.stats.dfe=grad_dfe;
stresults{i}.stats.yr=residuals(grad_lme);
stresults{i}.tstat=grad_t;
WM_regional_grad_t(i,1)=stresults{i}.tstat;
WM_regional_grad_t(i,2)=stresults{i}.pval;
end
save('G2_related_2_back_dprime_stresults','stresults');
save('G2_related_2_back_dprime_t','WM_regional_grad_t');


% corr WM and ROI with significant age effect on primary-transmodal gradient 

maskfile='...\grey_mask_0.2_4mm_AAL90_mask.nii';
mask_hdr = spm_vol(maskfile);
mask_vol = spm_read_vols(mask_hdr);
mask_ind = reshape(mask_vol>0,1,[]);
for i=1:6
ROIfile=['G1_ROI',num2str(i),'_mask.nii'];
ROI_hdr = spm_vol(ROIfile);
ROI_vol = spm_read_vols(ROI_hdr);
ROI_ind=ROI_vol(mask_ind);
index=logical(ROI_ind);
ROI_grad=indi_grad1(:,index);
ROI_grad1(:,i)=mean(ROI_grad,2);
end

ROI_grad1_wm=ROI_grad1(index,:);
for i=1:6
ROI=ROI_grad1_wm(:,i);
table_model_resort_WM.ROI=ROI; 
[ROI_t,ROI_pValue,ROI_beta,ROI_lme,ROI_dfe] = mixed_model_WM(ROI,table_model_resort_WM);
corr_ROI_grad1_WM(i,1)=ROI_t;
corr_ROI_grad1_WM(i,2)=ROI_pValue;
end

