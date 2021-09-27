%%examine the age-associated changes

% calculate the global gradient topography

load('realigned_gradient.mat')
indi_grad1(:,:)=realigned_gradient.realigned(:,1,:)';
indi_grad2(:,:)=realigned_gradient.realigned(:,2,:)';
range_indi_grad(:,1)=range(indi_grad1,2);
range_indi_grad(:,2)=range(indi_grad2,2);
std_indi_grad(:,1)=std(indi_grad1,0,2);
std_indi_grad(:,2)=std(indi_grad2,0,2);
save('indi_grad1','indi_grad1');
save('indi_grad2','indi_grad2');
save('range_indi_grad','range_indi_grad');
save('std_indi_grad','std_indi_grad');

for n=1:491
gradient(:,1)=indi_grad1(n,:)';
gradient(:,2)=indi_grad2(n,:)';
P(:,1)=median(gradient);
for j=1:length(gradient)
    P(:,2)=gradient(j,:)';
    Gdist_node=dist(P);
    Gdist(n,j)=Gdist_node(1,2);
end
end
save Gdist Gdist

for n=1:491
Gdist_global(n,1)=sumsqr(Gdist(n,:));
end
save Dispersion_Global Gdist_global

% aliged explanation ratio
gradient_order_resort=zeros(491,11); % 491: number of subjects
for i=1:491
    A=realigned_gradient.xfms{i};
    gradient_order=abs(A);
    [max_gradient_order,index]= max(gradient_order);
    gradient_order_resort(i,:)=index; 
end

gradient_indi_explan=zeros(491,136);
for i=1:491
load(['...\g',sublist{i}]);
explan=gradient.res.lambdas./sum(gradient.res.lambdas);
gradient_indi_explan(i,:)=explan;
end

gradient_indi_explan_resort=zeros(491,2); 
for i=1:491
gradient_indi_explan_resort(i,1)=gradient_indi_explan(i,gradient_order_resort(i,1));
gradient_indi_explan_resort(i,2)=gradient_indi_explan(i,gradient_order_resort(i,2));
end

save('gradient_order_resort','gradient_order_resort');
save('gradient_indi_explan','gradient_indi_explan');
save('gradient_indi_explan_resort','gradient_indi_explan_resort');

%mixed linear model

load('...\age_associated_ana\table_model_resort.mat');
for i=1:17673 
prediction=indi_grad1(:,i);
table_model_resort.depen_var=prediction;
lme1 = fitlme(table_model_resort,'depen_var ~ 1 + age + sex + meanFD + (1|subname) + (-1 + age|subname) ');  
lme2 = fitlme(table_model_resort,'depen_var ~ 1 + age^2 + sex + meanFD + (1|subname) + (-1 + age|subname) + (-1 - age + age^2|subname) ');  
if lme2.Coefficients.pValue(5)<0.05
    %AIC selecte age model
    if  lme2.ModelCriterion.AIC > lme1.ModelCriterion.AIC
        model_type = 1;
        age_pValue = lme1.Coefficients.pValue(2);
        age_beta = lme1.Coefficients.Estimate(2);
        age_t=lme1.Coefficients.tStat(2);
        age_lme=lme1;
        age_dfe=lme1.DFE;
    else
        model_type = 2;
        age_pValue = lme2.Coefficients.pValue(5);    %significance age^2
        age_beta = lme2.Coefficients.Estimate(5);
        age_t=lme2.Coefficients.tStat(5);
        age_lme=lme2;
        age_dfe=lme2.DFE;
    end
else
    model_type = 1;
    age_pValue = lme1.Coefficients.pValue(2);
    age_beta = lme1.Coefficients.Estimate(2);
    age_t=lme1.Coefficients.tStat(2);
    age_lme=lme1;
    age_dfe=lme1.DFE;
end
stresults{i}.b=age_beta;  
stresults{i}.pval=age_pValue;
stresults{i}.stats.dfe=age_dfe;
stresults{i}.stats.yr=residuals(age_lme);
t_test(i,1)=age_t;
t_test(i,2)=age_pValue;
end
save('age_related_indi_grad1_stresults','stresults');
save('age_related_indi_grad1_t','t_test');

% GRF correction
hdr_mask = spm_vol('...\grey_mask_0.2_4mm_AAL90_mask.nii');
vol_mask = spm_read_vols(hdr_mask);
ind = find(vol_mask);
hdr = hdr_mask;
hdr.dt(1) = 16;
vol = zeros(hdr.dim);
R_Volume = zeros([size(vol),491]); % subjects num 
DOF = 0;
for i = 1:length(ind)
    [x,y,z] = ind2sub(size(vol),ind(i));
    R_Volume(x,y,z,:) = stresults{i}.stats.yr;
    DOF = DOF + stresults{i}.stats.dfe;
end

MaskData = vol_mask;

DOF = DOF/length(ind);

Vox = [4 4 4];

R_Volume=single(R_Volume);
[N1, N2, N3, N4]=size(R_Volume);

R_Volume=(R_Volume-repmat(mean(R_Volume,4),[1,1,1, N4]))./repmat(std(R_Volume,0,4),[1,1,1, N4]);%Zero mean and one std
R_Volume(isnan(R_Volume))=0;

SSminus=[0 0 0];
S2=[0 0 0];

N=0;
for x=2:N1
    for y=2:N2
        for z=2:N3
            if MaskData(x, y, z) && MaskData(x-1, y, z) && MaskData(x, y-1, z) && MaskData(x, y, z-1)
                N=N+1;
                for t=1:N4
                    SSminus(1) = SSminus(1) + R_Volume(x, y, z, t) * R_Volume(x-1, y, z, t);
                    SSminus(2) = SSminus(2) + R_Volume(x, y, z, t) * R_Volume(x, y-1, z, t);
                    SSminus(3) = SSminus(3) + R_Volume(x, y, z, t) * R_Volume(x, y, z-1, t);
                    
                    S2(1) = S2(1) + 0.5 * ((R_Volume(x, y, z, t)^2) + (R_Volume(x-1, y, z, t)^2));
                    S2(2) = S2(2) + 0.5 * ((R_Volume(x, y, z, t)^2) + (R_Volume(x, y-1, z, t)^2));
                    S2(3) = S2(3) + 0.5 * ((R_Volume(x, y, z, t)^2) + (R_Volume(x, y, z-1, t)^2));
                end
            end
        end
    end
    %     fprintf('.');
end

if SSminus(1)>0.99999999*S2(1)
    SSminus(1)=0.99999999*S2(1);
    warning('possibly biased smootheness in X');
end
if SSminus(2)>0.99999999*S2(2)
    SSminus(2)=0.99999999*S2(2);
    warning('possibly biased smootheness in Y');
end
if SSminus(3)>0.99999999*S2(3)
    SSminus(3)=0.99999999*S2(3);
    warning('possibly biased smootheness in Z');
end

sigmasq(1) = -1 / (4 * log(abs(SSminus(1)/S2(1))));
sigmasq(2) = -1 / (4 * log(abs(SSminus(2)/S2(2))));
sigmasq(3) = -1 / (4 * log(abs(SSminus(3)/S2(3))));

dLh=((sigmasq(1)*sigmasq(2)*sigmasq(3))^-0.5)*(8^-0.5);

if N4 > 1
    fprintf('DLH %f voxels^-3 before correcting for temporal DOF\n',dLh);
    
    lut(6)   = 1.5423138; lut(7)   = 1.3757105; lut(8)   = 1.2842680;
    lut(9)   = 1.2272151; lut(10)  = 1.1885232; lut(11)  = 1.1606988;
    lut(12)  = 1.1398000; lut(13)  = 1.1235677; lut(14)  = 1.1106196;
    lut(15)  = 1.1000651; lut(16)  = 1.0913060; lut(17)  = 1.0839261;
    lut(18)  = 1.0776276; lut(19)  = 1.0721920; lut(20)  = 1.0674553;
    lut(21)  = 1.0632924; lut(26)  = 1.0483053; lut(31)  = 1.0390117;
    lut(41)  = 1.0281339; lut(51)  = 1.0219834; lut(61)  = 1.0180339;
    lut(71)  = 1.0152850; lut(81)  = 1.0132621; lut(91)  = 1.0117115;
    lut(101) = 1.0104851; lut(151) = 1.0068808; lut(201) = 1.0051200;
    lut(301) = 1.0033865; lut(501) = 1.0020191;
    
    y = lut(lut~=0);
    x = find(lut~=0);
    xi=[1:501];
    lut_interpolated=interp1(x,y,xi,'linear');
    
    if (DOF < 6)
        dLh=dLh * 1.1;
    elseif (DOF>500)
        dLh=dLh * (1.0321/DOF +1)^0.5;
    else
        retval=(lut_interpolated(floor(DOF)+1)-lut_interpolated(floor(DOF)))*(floor(DOF)+1-floor(DOF)) + ...
            lut_interpolated(floor(DOF)+1);
        dLh=dLh * retval^0.5;
    end
    
end

FWHM(1) =  sqrt(8 * log(2) * sigmasq(1));
FWHM(2) =  sqrt(8 * log(2) * sigmasq(2));
FWHM(3) =  sqrt(8 * log(2) * sigmasq(3));

resels = FWHM(1)*FWHM(2)*FWHM(3);
fprintf('\nFWHMx = %f voxels\nFWHMy = %f voxels\nFWHMz = %f voxels\n',FWHM(1),FWHM(2),FWHM(3));
FWHM=FWHM.*Vox;
fprintf('FWHMx = %f mm\nFWHMy = %f mm\nFWHMz = %f mm\n',FWHM(1),FWHM(2),FWHM(3));
nVoxels=length(find(MaskData));
fprintf('DLH = %f\nVOLUME = %d\nRESELS = %f\n',dLh,nVoxels,resels);

VoxelPThreshold = 0.001;
ClusterPThreshold = 0.05/2;

zThrd=norminv(1 - VoxelPThreshold/2);
D=3;
Em = nVoxels * (2*pi)^(-(D+1)/2) * dLh * (zThrd*zThrd-1)^((D-1)/2) * exp(-zThrd*zThrd/2);
EN = nVoxels * (1-normcdf(zThrd)); %In Friston et al., 1994, EN = S*Phi(-u). (Expectation of N voxels)  % K. Friston, K. Worsley, R. Frackowiak, J. Mazziotta, and A. Evans. Assessing the significance of focal activations using their spatial extent. Human Brain Mapping, 1:214?220, 1994.
Beta = ((gamma(D/2+1)*Em)/(EN)) ^ (2/D); % K. Friston, K. Worsley, R. Frackowiak, J. Mazziotta, and A. Evans. Assessing the significance of focal activations using their spatial extent. Human Brain Mapping, 1:214?220, 1994.

% Get the minimum cluster size
pTemp=1;
ClusterSize=0;
while pTemp >= ClusterPThreshold
    ClusterSize=ClusterSize+1;
    pTemp = 1 - exp(-Em * exp(-Beta * ClusterSize^(2/D))); %K. Friston, K. Worsley, R. Frackowiak, J. Mazziotta, and A. Evans. Assessing the significance of focal activations using their spatial extent. Human Brain Mapping, 1:214?220, 1994.
end

fprintf('%f voxels\n',ClusterSize);








