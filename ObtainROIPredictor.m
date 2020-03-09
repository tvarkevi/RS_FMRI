function [predictor_vector_right_ROI, predictor_vector_left_ROI] = ObtainROIPredictor(fileInfo, region_of_interest, all_subject_scans, varargin)

% Obtain predictor variables for left and right side of each seed-region
% for analyses
%
% Input arguments(region_of_interest: Specify name of the
% region-of-interest, e.g.:
%       region_of_interest = 'rAmygdala_LB.nii';
% Input arguments(all_subject_scans): All scans for this subject, inherited
% from VoxelWiseCorrelationMatrix
% OPTIONAL => input argument(varargin): Specify name of nuisance region of
% no interest, e.g.:
%       nuisance_region = 'rAmygdala_CM.nii'
% Subfunctions: -


base_nargin = 3;

% ----- Read probabilistic map of seed region-of-interest ----- %
H_ROI = spm_vol([fileInfo.base_dir '\' region_of_interest]);
[D_ROI, XYZ_ROI] = spm_read_vols(H_ROI);

% ----- Define each hemisphere seperately for the ROI ----- %
right_ROI_map = D_ROI;         left_ROI_map = D_ROI;

all_ROI_indices = find(D_ROI > 0);                                          
all_ROI_XYZ_of_indices = XYZ_ROI(:, all_ROI_indices);                       

right_indices_of_XYZ = find(all_ROI_XYZ_of_indices(1,:) > 0)';              
left_indices_of_XYZ = find(all_ROI_XYZ_of_indices(1,:) < 0)';               

right_ROI_indices = all_ROI_indices(right_indices_of_XYZ);                  
left_ROI_indices = all_ROI_indices(left_indices_of_XYZ);                    

right_ROI_map(left_ROI_indices) = 0;                                        
left_ROI_map(right_ROI_indices) = 0;                                        

% ----- Optional: Control for partial volume sharing (method 1) ----- %
if nargin > base_nargin
    % ----- Obtain probability maps of nuisance region-of-interest ----- %
    [right_NUI_map, left_NUI_map] = inner_partial_volume_correction(fileInfo, varargin{1});
    
    % ----- Obtain partial volume corrected non-zero indices ----- %
    right_ROI_non_zero_indices = find(right_ROI_map > right_NUI_map);
    left_ROI_non_zero_indices = find(left_ROI_map > left_NUI_map);
else
    % ----- Find indices of non-zero (probability) values -----%
    right_ROI_non_zero_indices = find(right_ROI_map > 0);
    left_ROI_non_zero_indices = find(left_ROI_map > 0);
end

% ----- Obtain predictor vector for the region-of-interest ----- %
predictor_vector_right_ROI = [];                                            
predictor_vector_left_ROI = [];                                             
    % ----- Obtain ROI weightsum for each scan of this subject ----- %
for iScan = 1:length(all_subject_scans)                        
    this_subjects_scan = all_subject_scans{iScan};                          
    
    right_hemi_data = this_subjects_scan;                                    
    left_hemi_data = this_subjects_scan;                                    
    
    % ----- Multiply BOLD with probability for this scan ----- %
    right_hemi_data(right_ROI_non_zero_indices) = right_hemi_data(right_ROI_non_zero_indices) .* right_ROI_map(right_ROI_non_zero_indices);      
    left_hemi_data(left_ROI_non_zero_indices) = left_hemi_data(left_ROI_non_zero_indices) .* left_ROI_map(left_ROI_non_zero_indices);           
    
    % ----- Obtain weighted sum of probability-weighted BOLD values ----- %
    right_ROI_weightsum = sum(right_hemi_data(right_ROI_non_zero_indices)); 
    left_ROI_weightsum = sum(left_hemi_data(left_ROI_non_zero_indices));    
    
    % ----- Add weightsum to the ROI predictor vector ----- %
    predictor_vector_right_ROI = [predictor_vector_right_ROI; right_ROI_weightsum]; 
    predictor_vector_left_ROI = [predictor_vector_left_ROI; left_ROI_weightsum];    
end

end

function [right_NUI_map, left_NUI_map] = inner_partial_volume_correction(fileInfo, nuisance_region)

% ----- Read probabilistic map of nuisance region ----- %
H_NUI = spm_vol([fileInfo.base_dir '\' nuisance_region]);
[D_NUI, XYZ_NUI] = spm_read_vols(H_NUI);

% ----- Define each hemisphere seperately for the nuisance variable ----- %
right_NUI_map = D_NUI;         left_NUI_map = D_NUI;

all_NUI_indices = find(D_NUI > 0);                                          
all_NUI_XYZ_of_indices = XYZ_NUI(:, all_NUI_indices);                       

right_indices_of_XYZ = find(all_NUI_XYZ_of_indices(1,:) > 0)';              
left_indices_of_XYZ = find(all_NUI_XYZ_of_indices(1,:) < 0)';               

right_NUI_indices = all_NUI_indices(right_indices_of_XYZ);                  
left_NUI_indices = all_NUI_indices(left_indices_of_XYZ);                    

right_NUI_map(left_NUI_indices) = 0;                                        
left_NUI_map(right_NUI_indices) = 0;

end
