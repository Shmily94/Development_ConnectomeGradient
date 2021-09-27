% calculate long- and short-range FCS
for n=1:491
load(sublist(n).name);
thresholdMatrix = threshold_proportional(zFC,0.1);
Edist_thresh_ind=logical(thresholdMatrix);
A=double(Edist_thresh_ind);
Edist_thresh=Edist.* A; %

Dist_higher_70_ind=find(Edist_thresh>70);
FC_dist_higher_70(n,2)=mean(Edist_thresh(Dist_higher_70_ind));
FC_dist_higher_70(n,1)=mean(thresholdMatrix(Dist_higher_70_ind));

Dist_lower_70_ind=find(Edist_thresh>0 & Edist_thresh<70);
FC_dist_lower_70(n,2)=mean(Edist_thresh(Dist_lower_70_ind));
FC_dist_lower_70(n,1)=mean(thresholdMatrix(Dist_lower_70_ind));

Dist_70_ind=find(Edist_thresh>60 & Edist_thresh<70);
FC_dist_70(n,2)=mean(Edist_thresh(Dist_70_ind));
FC_dist_70(n,1)=mean(thresholdMatrix(Dist_70_ind));

Dist_60_ind=find(Edist_thresh>50 & Edist_thresh<60);
FC_dist_60(n,2)=mean(Edist_thresh(Dist_60_ind));
FC_dist_60(n,1)=mean(thresholdMatrix(Dist_60_ind));

Dist_50_ind=find(Edist_thresh>40 & Edist_thresh<50);
FC_dist_50(n,2)=mean(Edist_thresh(Dist_50_ind));
FC_dist_50(n,1)=mean(thresholdMatrix(Dist_50_ind));

Dist_40_ind=find(Edist_thresh>30 & Edist_thresh<40);
FC_dist_40(n,2)=mean(Edist_thresh(Dist_40_ind));
FC_dist_40(n,1)=mean(thresholdMatrix(Dist_40_ind));

Dist_30_ind=find(Edist_thresh>20 & Edist_thresh<30);
FC_dist_30(n,2)=mean(Edist_thresh(Dist_30_ind));
FC_dist_30(n,1)=mean(thresholdMatrix(Dist_30_ind));

for k=1:8
    index=find(Network_8==k);
    A_index=A(index,index); 
    kwithin(n,k)=sum(A_index,'all');
    between_index=find(Network_8>0 & Network_8~=k);
    B_index=A(index,between_index);
    kbetween(n,k)=sum(B_index,'all');
    MSI(n,k)=(kwithin(n,k)-kbetween(n,k))/kwithin(n,k);
end
   
clear zFC thresholdMatrix Edist_thresh_ind 
end

save('MSI','MSI');


%corr network topology and global gradient topolography
load('...\table_model_resort');
table_model_resort.cov_age=table_model_resort.age;
table_model_resort.age=table_model_resort.network;
for i=1:8
network=MSI(:,i);
for j=1:3
prediction=global_gradient(:,j);
table_model_resort.depen_var=prediction;  
table_model_resort.network=network;
lme1 = fitlme(table_model_resort,'depen_var ~ 1 + network + sex + meanFD + cov_age + (1|subname) + (-1 + network|subname) ');  
lme2 = fitlme(table_model_resort,'depen_var ~ 1 + network^2 + sex + meanFD + cov_age + (1|subname) + (-1 + network|subname) + (-1 - network + network^2|subname) ');  
if lme2.Coefficients.pValue(5)<0.05
    %AIC selecte age model
    if  lme2.ModelCriterion.AIC > lme1.ModelCriterion.AIC
        model_type = 1;
        network_pValue = lme1.Coefficients.pValue(2);
        network_beta = lme1.Coefficients.Estimate(2);
        network_t=lme1.Coefficients.tStat(2);
        network_lme=lme1;
        network_dfe=lme1.DFE;
    else
        model_type = 2;
        network_pValue = lme2.Coefficients.pValue(5);    %significance age^2
        network_beta = lme2.Coefficients.Estimate(5);
        network_t=lme2.Coefficients.tStat(5);
        network_lme=lme2;
        network_dfe=lme2.DFE;
    end
else
    model_type = 1;
    network_pValue = lme1.Coefficients.pValue(2);
    network_beta = lme1.Coefficients.Estimate(2);
    network_t=lme1.Coefficients.tStat(2);
    network_lme=lme1;
    network_dfe=lme1.DFE;
end
stresults_MSI__grad{i,j}.b=network_beta;  
stresults_MSI__grad{i,j}.pval=network_pValue;
stresults_MSI__grad{i,j}.stats.dfe=network_dfe;
stresults_MSI__grad{i,j}.stats.yr=residuals(network_lme);
stresults_MSI__grad{i,j}.tstat=network_t;
end
end
save('MSI_related_grad','stresults_MSI__grad');
