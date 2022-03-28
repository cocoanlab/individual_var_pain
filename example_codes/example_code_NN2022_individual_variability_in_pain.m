%% Example code for "Individual Variaibility in Brain Representations of Pain", 
%  Kohoutova et al., 2022, Nature Neuroscience

% This code provides an example workflow of the multivariate analysis presented in the paper.
% The steps below include:
% 1. Building individualised models + permutation test
% 2. Apply masks to individualised models 
% 3. Calculate correlation between individualised regional patterns 
% 4. Get permuted correlation distance
% 5. Calculate z-scores (normalised representational distance)

%% Set directories and get file names and region masks
datdir = 'betamaps/directory';
savedir_model = 'directory/to/save/individualised/models';
savedir_perm = 'directory/to/save/permuted/data';


allfiles = dir(fullfile(datdir, 'sub*.mat')); % subject data in fmri_data format
allfiles = {allfiles.name};

maskdir = 'directory/containing/region/masks';

masks = filenames(fullfile(maskdir, '*.nii'));

% load masks
all_masks = cell(21,1);
for m = 1:length(masks.fnames)

            temp_mask = fmri_data(masks{m}, 'brainmask.nii');
            all_masks{m} = temp_mask.dat;

end

%% 1. Building individualised models + permutation test

for i = 1:length(allfiles) 
    
    fprintf('\nWorking on subject %02d/%02d', i, length(allfiles));
    
    % load subject's data
    load(fullfile(datdir, allfiles{i})) 
    
    % train individualised SVR model with 5-fold CV
    [~, out] = predict(sub_dat, 'algorithm_name', 'cv_svr', 'nfolds', 5, 'error_type', 'mse');
    
    % for boostrap with 10,000 samples, run 
    % [~, boot] = predict(sub_dat, 'algorithm_name', 'cv_svr', 'nfolds', 5, 'error_type', 'mse', 'bootsamples', 5000);
    
    sub_id = split(allfiles{i}, ".");
    sub_id = sub_id{1};
    
    savename = fullfile(savedir_model, sub_id);
    save(savename, 'out')
    
    % run permutation - 1,000 iterations; save only regional data 
    perm_subj_masked = cell(21,1);
    fprintf('Permutation start \n');
    for j = 1:1000
        
        % randomise ratings
        rand_y_id = randperm(length(sub_dat.Y));
        rand_y = sub_dat.Y(rand_y_id);

        dat_temp = sub_dat;
        dat_temp.Y = rand_y;

        % train SVR model; no CV
        [~, perm] = predict(dat_temp, 'algorithm_name', 'cv_svr', 'nfolds', 1, 'error_type', 'mse');
        
        % get weight map only
        map = perm.weight_obj;
        
        % normalise weight map
        norm_dat = map.dat./std(map.dat);
        
        % divide into regions
        for k = 1:length(all_masks)
                
                perm_temp = norm_dat.*all_masks{k};
                perm_subj_masked{k}(:,j) = perm_temp(all_masks{k}==1);
                
        end

    end
    
    savename = fullfile(savedir_perm, filename);
    save(savename, 'perm_subj_masked')

    
end

%% 2. Apply masks to individualised models

modelnames= filenames(fullfile(savedir_model, '*.mat'));

% save all maps in one matrix
all_dat = zeros(328798,124);
for i = 1:length(modelnames)
    
    load(modelnames{i})
    dat_temp = out.weight_obj.dat;
    all_dat(:,i) = dat_temp;
    
end

% normalise
norm_dat = all_dat ./ repmat(std(all_dat), size(all_dat,1), 1);

dat_temp_masked = zeros(size(norm_dat,1),size(norm_dat,2));
dat_all_masked = cell(length(masks),1);

% mask data
for i = 1:length(masks)
    
    temp_mask = all_masks{i};
    
        for ii = 1:size(norm_dat,2) % for each subject
            
        dat_temp_masked(:,ii) = norm_dat(:,ii).*temp_mask;
        
        end
        
    dat_all_masked{i} = dat_temp_masked;
end

%% 3. Calculate correlation between individualised regional patterns 

% number of subjects
n_subject = 404; 

sim_mat = zeros(nchoosek(n_subject,2),1);
corr_mat_all = zeros(nchoosek(n_subject,2),length(masks));

for j = 1:length(dat_all_masked)
    
     
    fprintf('\nWorking on region %02d/%02d', j, length(masks));
    k =1;

    dat_masked_temp = dat_all_masked{j};
    for i = 1:n_subject
        
        
        for ii = 1:n_subject
            
            if i<ii
            sim_mat(k) = corr(dat_masked_temp(:,i),dat_masked_temp(:,ii));
            k = k+1;
            end
        end
    end
    
    corr_mat_all(:,j) = sim_mat;

end

%% 4. Get permuted correlation distance

all_perm_dat = filenames(fullfile(savedir_perm, '*.mat'));

perm_corr_dist_all = cell(length(masks),1);
m = 1;
for i = 1:n_subject 
    
    fprintf('\nWorking on subject %02d/%02d', i, n_subject);
    
    load(all_perm_dat{i})
    subj1 = perm_subj_masked;
    
    for j = 1:n_subject 
        
        if i<j 
            
            fprintf('\nwith subject %02d/%02d', j, n_subject);
            load(all_perm_dat{j})
            subj2 = perm_subj_masked;
            
            for k = 1:length(masks)
                
                fprintf('\nand region %02d/%02d', k, length(masks));
                reg_temp1 = double(subj1{k});
                reg_temp2 = double(subj2{k});
                
                reg_corr_temp = 1-col_corr(reg_temp1, reg_temp2); % column correlation defined below
                
                perm_corr_dist_all{k}(m,:) = reg_corr_temp;
               
                
            end
            m = m+1;
        end
        
    end

end

%% 5. Calculate z-scores (normalised representational distance)

% observed correlation to correlation distance (representational dissimilarity)
corr_mat_all_dist = 1-corr_mat_all; 

zscore_corr_dist = zeros(nchoosek(n_subject,2),length(masks));
for i = 1:21
    
    fprintf('\nWorking on region %02d/%02d', i, length(masks));
    
    % permuted dissimilarity matrix - null-hypothesis baseline 
    perm_dat_temp = perm_corr_dist_all{i};
    
    std_corrd = std(perm_dat_temp,[],2);
    mean_corrd = mean(perm_dat_temp,2);
    
    % calculate z-scores
    zscore_corr_dist(:,i) = (corr_mat_all_dist(:,i)-mean_corrd)./std_corrd;
    
    
    
end





%% Column correlation function

function C = col_corr(A, B)


A = bsxfun(@minus, A, mean(A));
B = bsxfun(@minus, B, mean(B));

A = bsxfun(@times, A, 1./sqrt(sum(A.^2)));
B = bsxfun(@times, B, 1./sqrt(sum(B.^2)));

C = sum(A.*B);

end
