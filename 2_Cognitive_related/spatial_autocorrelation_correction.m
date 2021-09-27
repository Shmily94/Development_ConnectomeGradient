%% generate surrogate t-maps for corgnitive terms
hdr_mask = spm_vol('...\grey_mask_0.2_4mm_AAL90_mask.nii');
vol_mask = spm_read_vols(hdr_mask);
ind = find(vol_mask);
nVoxels = length(ind);
real_g1_t_hdr = spm_vol('g1_tstat.nii');

tok=regexp(real_g1_t_hdr.descrip, '\{dLh_(.*?)\}\{FWHMx_(.*?)FWHMy_(.*?)FWHMz_(.*?)mm\}',...
    'tokens');
if isempty(tok) || numel(tok)~=1
    return;
end
dLh=str2double(tok{1}{1});

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


load('surrogate_maps_g1_t_resample.mat');

parfor j = 1:10000
    id = ['0000',num2str(j)];
    pos_flag = 0;
    neg_flag = 0;
    cd ...\neurosynth_analysis\surr_g1_pos
    if exist(['rsurr_g1_t_pos_',id(end-4:end),'.nii'],'file')
        pos_flag = 1;
    end
    cd ...\neurosynth_analysis\surr_g1_neg
    if exist(['rsurr_g1_t_neg_',id(end-4:end),'.nii'],'file')
        neg_flag = 1;
    end
        
    if pos_flag*neg_flag == 0
        vol_surr = zeros(hdr_mask.dim);
        vol_surr(ind) = surrogate_maps(j,:);
        vol_surr(vol_surr < zThrd & vol_surr > -zThrd) = 0;
        [L, num] = bwlabeln(vol_surr,26);
        n = 0;
        for x = 1:num
            theCurrentCluster = L == x;
            if length(find(theCurrentCluster)) <= ClusterSize
                n = n + 1;
                vol_surr(logical(theCurrentCluster)) = 0;
            end
        end
        if pos_flag == 0
            vol_pos = vol_surr;
            vol_pos(vol_surr<0) = 0;
            hdr_pos = real_g1_t_hdr;
            hdr_pos.fname = ['surr_g1_t_pos_',id(end-4:end),'.nii'];
            cd ...\neurosynth_analysis\surr_g1_pos
            spm_write_vol(hdr_pos,vol_pos);
            x_reslice('...\neurosynth_analysis\rG1_t_pos.nii',['surr_g1_t_pos_',id(end-4:end),'.nii'],1);
            delete(['surr_g1_t_pos_',id(end-4:end),'.nii']);
        end
        
        if neg_flag == 0
            vol_neg = -vol_surr;
            vol_neg(vol_surr>0) = 0;
            hdr_neg = real_g1_z_hdr;
            hdr_neg.fname = ['surr_g1_t_neg_',id(end-4:end),'.nii'];
            cd ...\neurosynth_analysis\surr_g1_neg
            spm_write_vol(hdr_neg,vol_neg);
            x_reslice('...\neurosynth_analysis\rG1_t_pos.nii',['surr_g1_t_neg_',id(end-4:end),'.nii'],1);
            delete(['surr_g1_t_neg_',id(end-4:end),'.nii']);
        end
    end
end
% spatial autocorrelation correction
real_r = zeros(30,1);
fid = fopen('real_r_g1_pos.txt');
c = textscan(fid, '%s %s','delimiter',',');
fclose(fid);
real_r = str2double(c{1,2}(2:end));

surr_r = zeros(30,10000);
cd ...\neurosynth_analysis\surr_g1_pos
parfor i = 1:10000
    id = ['0000',num2str(i)];
    fid = fopen(['rsurr_g1_t_pos_',id(end-4:end),'.txt']);
    c = textscan(fid, '%s %s','delimiter',',');
    fclose(fid);
    surr_r(:,i) = str2double(c{1,2}(2:end));
end

p_cogn_term = sum(gt(surr_r,real_r),2)/10000;
cd ...\neurosynth_analysis
save cog_term_g1_pos.mat real_r surr_r p_cogn_term

real_r = zeros(30,1);
fid = fopen('real_r_g1_neg.txt');
c = textscan(fid, '%s %s','delimiter',',');
fclose(fid);
real_r = str2double(c{1,2}(2:end));

surr_r = zeros(30,10000);
cd ...\neurosynth_analysis\surr_g1_neg
parfor i = 1:10000
    id = ['0000',num2str(i)];
    fid = fopen(['rsurr_g1_t_neg_',id(end-4:end),'.txt']);
    c = textscan(fid, '%s %s','delimiter',',');
    fclose(fid);
    surr_r(:,i) = str2double(c{1,2}(2:end));
end

p_cogn_term = sum(gt(surr_r,real_r),2)/10000;
cd ...\neurosynth_analysis
save cog_term_g1_neg.mat real_r surr_r p_cogn_term
