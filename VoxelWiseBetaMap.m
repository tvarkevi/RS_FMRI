function [correlation_map_right_ROI, correlation_map_left_ROI, H_Scan] = VoxelWiseBetaMap(seed_region, iSubject, fileInfo, CSF_parameters, white_matter_parameters, varargin)

% Conduct regression analyses (GLM): correlate the seed-region to each
% voxel in the brain (i.e., voxel-wize regression).
%
% input argument(seed_region): specify seed-region to base the connectivity
% analysis on, e.g.,
%       seed_region = 'rAmygdala_LB.nii'
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


base_nargin = 5;

% ----- Read functional scans of this subject ----- %
nScans = length(fileInfo.functional_file_names{iSubject});
for iScan = 1:nScans
    H_Scan = spm_vol(fileInfo.functional_file_names{iSubject}{iScan});
    D_Scan = spm_read_vols(H_Scan);
    all_subject_scans{iScan} = D_Scan;
end 

% ----- Obtain predictor vector of seed region ----- %
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

% ----- Compute empty correlation matrices ----- %
correlation_map_right_ROI = zeros(size(all_subject_scans{1}));                                                  
correlation_map_left_ROI = zeros(size(all_subject_scans{1})); 

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

% ----- Compute time-courses of voxels ----- %
for iVoxel = 1:length(find_brain_voxels)                                                                                        
    this_voxel_index = find_brain_voxels(iVoxel);
    
    [x, y, z] = ind2sub(size(all_subject_scans{1}), this_voxel_index);                         
    this_voxel_signal = [];                                                 
    for iScan = 1:length(all_subject_scans)                                             
        this_voxel_signal = [this_voxel_signal; all_subject_scans{iScan}(x, y, z)];        
    end
    
    % ----- Filter voxel signal ----- %
    [~, this_voxel_signal] = ezfilt(this_voxel_signal, 1/TR, 0.01);           
    [this_voxel_signal, ~] = ezfilt(this_voxel_signal, 1/TR, 0.08);           
    
    % ----- Specify regression model ----- %
    outcome_variable = this_voxel_signal(:);                                                    
    outcome_variable = zscore(outcome_variable);                                                              
    
    if nargin > base_nargin
        predictor_variables = [predictor_vector_right_ROI(:) predictor_vector_right_NUI(:) global_mean(:) fileInfo.motion_parameters{iSubject} white_matter_parameters CSF_parameters];
    else
        predictor_variables = [predictor_vector_right_ROI(:) global_mean(:) fileInfo.motion_parameters{iSubject} white_matter_parameters CSF_parameters];
    end
    predictor_variables = zscore(predictor_variables);
    this_beta_estimate_right = inv(predictor_variables' * predictor_variables) * predictor_variables' * outcome_variable;                                                                             
    
    if nargin > base_nargin
        predictor_variables = [predictor_vector_left_ROI(:) predictor_vector_left_NUI(:) global_mean(:) fileInfo.motion_parameters{iSubject} white_matter_parameters CSF_parameters];
    else
        predictor_variables = [predictor_vector_left_ROI(:) global_mean(:) fileInfo.motion_parameters{iSubject} white_matter_parameters CSF_parameters];
    end
    predictor_variables = zscore(predictor_variables);
    this_beta_estimate_left = inv(predictor_variables' * predictor_variables) * predictor_variables' * outcome_variable;
    
    % ----- Fill empty correlation matrices ----- %
    predictor_of_interest_column = 1;
    correlation_map_right_ROI(this_voxel_index) = this_beta_estimate_right(predictor_of_interest_column);                             
    correlation_map_left_ROI(this_voxel_index) = this_beta_estimate_left(predictor_of_interest_column);
end
