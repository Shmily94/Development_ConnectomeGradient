clear,clc
load('.../100DS360scaledRobustSigmoidNSGRNAseqQC1LRcortex_ROI_NOdistCorrEuclidean.mat','parcelExpression')
load('.../100DS360scaledRobustSigmoidNSGRNAseqQC1LRcortex_ROI_NOdistCorrEuclidean.mat','probeInformation')
load('.../beta_glasser.mat')

temp1=find(isnan(beta_glasser_z));
temp2=find(isnan(parcelExpression(:,2)));
missingdata_regions=union(temp1,temp2);
region_ind=setdiff(parcelExpression(:,1),missingdata_regions);

group_express=parcelExpression(region_ind,2:end);
gene_name = probeInformation.GeneSymbol;


[r p]=corr(beta_glasser(region_ind,:),group_express);
[r_rank ind]=sort(r,'descend');
gene_name_rank=gene_name(ind);

output_dir= '...\gene_ana';

fid = fopen(fullfile(output_dir,'corr_geneWeights.csv'),'w');
for i=1:length(gene_name_rank)
  fprintf(fid,'%s, %d, %f\n', gene_name_rank{i}, r_rank(i),r_rank(i));
end
fclose(fid);
