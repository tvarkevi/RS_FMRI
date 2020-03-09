function betas_seed_vs_ROI = ROIWiseGLM(seed_region, non_seed_region, iSubject, fileInfo, CSF_parameters, white_matter_parameters, varargin)

% Conduct regression analyses (GLM): correlate the seed-region to the
% region-of-interest, as specified by the input argument 'non_seed_region'.
%
% input argument(seed_region): specify seed-region to base the connectivity
% analysis on, e.g.,
%       seed_region = 'rAmygdala_LB.nii'
% input argument(non_seed_region): specify non-seed region-of-interest to 
% base the connectivity analysis on. Note: The ROI map needs to contain
% both hemispheres! E.g.,
%       non_seed_region = 'rOFC_Total.nii'
% Input argument(functional_file_names): specify full directories for each
% scan in timeseries.
% Input argument(motion_parameters): input data matrix that contains the
% estimates for motion (as nuisance variable).
% Input argument(CSF_parameters): input data matrix that contains the
% estimates for CSF (as nuisance variable).
% Input argument(white_matter_parameters): input data matrix that contains the
% estimates for white matter (as nuisance variable).
% OPTIONAL => input argument(varargin): specify nuisance region of no
% interest, e.g.:
%       nuisance_region = 'rAmygdala_CM.nii'
% Subfunctions: ObtainROIPredictor, ezfilt


base_nargin = 6;

% ----- Read functional scans of this subject ----- %
nScans = length(fileInfo.functional_file_names{iSubject});
for iScan = 1:nScans
    H_Scan = spm_vol(fileInfo.functional_file_names{iSubject}{iScan});
    D_Scan = spm_read_vols(H_Scan);
    all_subject_scans{iScan} = D_Scan;
end 

% ----- Obtain predictor vector of seed ROI ----- %
if nargin > base_nargin
    % ----- Optional: Control for partial volume sharing (method 1) ----- %
    [predictor_vector_right_ROI, predictor_vector_left_ROI] = ObtainROIPredictor(fileInfo, seed_region, all_subject_scans, varargin{1});
    % ----- Optional: Control for partial volume sharing (method 2) ----- %
    [predictor_vector_right_NUI, predictor_vector_left_NUI] = ObtainROIPredictor(fileInfo, varargin{1}, all_subject_scans, seed_region);
else
    [predictor_vector_right_ROI, predictor_vector_left_ROI] = ObtainROIPredictor(fileInfo, seed_region, all_subject_scans);
end

% ----- Filter signals of predictor variables ----- %
TR = 1.6;

[~, predictor_vector_right_ROI] = ezfilt(predictor_vector_right_ROI, 1/TR, 0.01);
[predictor_vector_right_ROI, ~] = ezfilt(predictor_vector_right_ROI, 1/TR, 0.08);

[~, predictor_vector_left_ROI] = ezfilt(predictor_vector_left_ROI, 1/TR, 0.01);
[predictor_vector_left_ROI, ~] = ezfilt(predictor_vector_left_ROI, 1/TR, 0.08);

if nargin > base_nargin
    [~, predictor_vector_right_NUI] = ezfilt(predictor_vector_right_NUI, 1/TR, 0.01);
    [predictor_vector_right_NUI, ~] = ezfilt(predictor_vector_right_NUI, 1/TR, 0.08);
    
    [~, predictor_vector_left_NUI] = ezfilt(predictor_vector_left_NUI, 1/TR, 0.01);
    [predictor_vector_left_NUI, ~] = ezfilt(predictor_vector_left_NUI, 1/TR, 0.08);
end

[~, white_matter_parameters] = ezfilt(white_matter_parameters, 1/TR, 0.01);
[white_matter_parameters, ~] = ezfilt(white_matter_parameters, 1/TR, 0.08);

[~, CSF_parameters] = ezfilt(CSF_parameters, 1/TR, 0.01);
[CSF_parameters, ~] = ezfilt(CSF_parameters, 1/TR, 0.08);

% ----- Find voxels corresponding to brain-tissue ----- %
file_parts = strsplit(fileInfo.anatomical_file_names{iSubject}, '\');
subject_directory = strjoin(file_parts(1:end-1), '\');
H_Mask = spm_vol([subject_directory '\Brain_Mask.nii']);
D_Mask = spm_read_vols(H_Mask);
find_brain_voxels = find(D_Mask > 0);

% ----- Compute global-mean signal and filter ----- % 
global_mean = [];                                                           
for iScan = 1:length(all_subject_scans)                                                  
    global_mean = [global_mean; mean(all_subject_scans{iScan}(find_brain_voxels))];    
end
[~, global_mean] = ezfilt(global_mean, 1/TR, 0.01);                       
[global_mean, ~] = ezfilt(global_mean, 1/TR, 0.08);                     

% ----- Obtain outcome vector of (non-seed) ROI, and filter signal ----- %
[outcome_vector_right_ROI, outcome_vector_left_ROI] = ObtainROIPredictor(fileInfo, non_seed_region, all_subject_scans);

[~, outcome_vector_right_ROI] = ezfilt(outcome_vector_right_ROI, 1/TR, 0.01);
[outcome_vector_right_ROI, ~] = ezfilt(outcome_vector_right_ROI, 1/TR, 0.08);

[~, outcome_vector_left_ROI] = ezfilt(outcome_vector_left_ROI, 1/TR, 0.01);
[outcome_vector_left_ROI, ~] = ezfilt(outcome_vector_left_ROI, 1/TR, 0.08);

% ----- Specify regression models: Preperation ----- %
outcome_vector_right_ROI = zscore(outcome_vector_right_ROI);
outcome_vector_left_ROI = zscore(outcome_vector_left_ROI);

predictor_of_interest_column = 1;

% ----- Specify regression models: right (non-seed) ROI ----- %
if nargin > base_nargin
    predictor_variables = [predictor_vector_right_ROI(:) predictor_vector_right_NUI(:) global_mean(:) fileInfo.motion_parameters{iSubject} white_matter_parameters CSF_parameters];
else
    predictor_variables = [predictor_vector_right_ROI(:) global_mean(:) fileInfo.motion_parameters{iSubject} white_matter_parameters CSF_parameters];
end
predictor_variables = zscore(predictor_variables);
beta_estimates = inv(predictor_variables' * predictor_variables) * predictor_variables' * outcome_vector_right_ROI;
betas_seed_vs_ROI{1, 1} = 'Right_seed_vs_right_ROI';
betas_seed_vs_ROI{1, 2} = beta_estimates(predictor_of_interest_column);

if nargin > base_nargin
    predictor_variables = [predictor_vector_left_ROI(:) predictor_vector_left_NUI(:) global_mean(:) fileInfo.motion_parameters{iSubject} white_matter_parameters CSF_parameters];
else
    predictor_variables = [predictor_vector_left_ROI(:) global_mean(:) fileInfo.motion_parameters{iSubject} white_matter_parameters CSF_parameters];
end
predictor_variables = zscore(predictor_variables);
beta_estimates = inv(predictor_variables' * predictor_variables) * predictor_variables' * outcome_vector_right_ROI;
betas_seed_vs_ROI{2, 1} = 'Left_seed_vs_right_ROI';
betas_seed_vs_ROI{2, 2} = beta_estimates(predictor_of_interest_column);

% ----- Specify regression models: left (non-seed) ROI ----- %
if nargin > base_nargin
    predictor_variables = [predictor_vector_right_ROI(:) predictor_vector_right_NUI(:) global_mean(:) fileInfo.motion_parameters{iSubject} white_matter_parameters CSF_parameters];
else
    predictor_variables = [predictor_vector_right_ROI(:) global_mean(:) fileInfo.motion_parameters{iSubject} white_matter_parameters CSF_parameters];
end
predictor_variables = zscore(predictor_variables);
beta_estimates = inv(predictor_variables' * predictor_variables) * predictor_variables' * outcome_vector_left_ROI;
betas_seed_vs_ROI{3, 1} = 'Right_seed_vs_left_ROI';
betas_seed_vs_ROI{3, 2} = beta_estimates(predictor_of_interest_column);

if nargin > base_nargin
    predictor_variables = [predictor_vector_left_ROI(:) predictor_vector_left_NUI(:) global_mean(:) fileInfo.motion_parameters{iSubject} white_matter_parameters CSF_parameters];
else
    predictor_variables = [predictor_vector_left_ROI(:) global_mean(:) fileInfo.motion_parameters{iSubject} white_matter_parameters CSF_parameters];
end
predictor_variables = zscore(predictor_variables);
beta_estimates = inv(predictor_variables' * predictor_variables) * predictor_variables' * outcome_vector_left_ROI;
betas_seed_vs_ROI{4, 1} = 'Left_seed_vs_left_ROI';
betas_seed_vs_ROI{4, 2} = beta_estimates(predictor_of_interest_column);
