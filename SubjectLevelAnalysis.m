function SubjectLevelAnalysis(seed_region, non_seed_region, analysis_type, fileInfo, varargin)

% Conduct first (subject) level analysis: i.e., resting-state connectivity
% anslysis.
%
% Input argument(seed_region): Specify seed-region to base the connectivity
% analysis on, e.g.,
%       seed_region = 'rAmygdala_Total.nii';
% Input argument(non_seed_region): Specify non-seed region-of-interest to
% base the connectivity analysis on, e.g.,
%       non_seed_region = 'rOFC_total.nii';
% Input argument(analysis_type): Specify if voxel-wise or ROI-wise analysis
% is to be conducted
%       1. Voxel-wise (exploratory) analysis
%       2. ROI-wise (confirmatory) analysis
% Input argument(fileInfo)
% OPTIONAL => input argument(varargin): specify nuisance region of no
% interest, e.g.:
%       nuisance_region = 'rAmygdala_CM.nii'
% Subfunctions: 
%       ObtainWhiteMatterPredictor/ObtainCSFPredictor (Optional)
%       VoxelWiseCorrelationMatrix/ROIWiseCorrelation


base_nargin = 4;

% ----- Nuisance variables: obtain white matter and CSF signal ----- %
% white_matter_parameters = ObtainWhiteMatterPredictor(fileInfo, 0);
% CSF_parameters = ObtainCSFPredictor(fileInfo, 0);
load('white_matter_parameters.mat');
load('CSF_parameters.mat');

% ----- First-level analysis: Loop over subjects ----- %
for iSubject = 1:length(fileInfo.functional_file_names)
    [save_dir, ~, ~] = fileparts(fileInfo.functional_file_names{iSubject}{1});
    
    % ----- Display progress (this subject) ----- %
    disp([num2str(iSubject) 9 save_dir 9 9 datestr(now)]); % The ASCI character for tab is 9.    
    
    switch analysis_type
        case 1 % Voxel-wise (exploratory) anslysis
            % ----- Obtain correlation matrices for both hemispheres of ROI ----- %
            if nargin > base_nargin
                [correlation_matrix_right_ROI, correlation_matrix_left_ROI, H_Scan] = VoxelWiseBetaMap(seed_region, ...
                    iSubject, fileInfo, CSF_parameters{iSubject}, white_matter_parameters{iSubject}, varargin{1});
            else
                [correlation_matrix_right_ROI, correlation_matrix_left_ROI, H_Scan] = VoxelWiseBetaMap(seed_region, ...
                    iSubject, fileInfo, CSF_parameters{iSubject}, white_matter_parameters{iSubject});
            end
            
            % ----- Set voxels outside of brain to not-a-number ----- %
            H_Mask = spm_vol([fileparts(fileInfo.anatomical_file_names{iSubject}) '\Brain_Mask.nii']);
            D_Mask = spm_read_vols(H_Mask);
            non_brain_voxels = find(D_Mask < 0.5);
            correlation_matrix_left_ROI(non_brain_voxels) = NaN;
            correlation_matrix_right_ROI(non_brain_voxels) = NaN;
            
            % ----- Save correlation matrix to (.nii) file ----- %
            H_Scan.fname = [save_dir '\Beta_map_left_' seed_region(1:end-4) '.nii'];
            H_Scan.pinfo(1) = 1;
            H_Scan.dt = [16 0];
            spm_write_vol(H_Scan, correlation_matrix_left_ROI);
            
            H_Scan.fname = [save_dir '\Beta_map_right_' seed_region(1:end-4) '.nii'];
            spm_write_vol(H_Scan, correlation_matrix_right_ROI);
            
        case 2 % ROI (confirmatory) analysis
            % ----- Obtain all combinations of seed vs. ROI models ----- %
            if nargin > base_nargin
                % ----- Optional: control for partial voluming ----- %
                betas_seed_vs_ROI = ROIWiseGLM(seed_region, non_seed_region, ...
                    iSubject, fileInfo, CSF_parameters{iSubject}, white_matter_parameters{iSubject}, ...
                    varargin{1});
            else
                betas_seed_vs_ROI = ROIWiseGLM(seed_region, non_seed_region, ...
                    iSubject, fileInfo, CSF_parameters{iSubject}, white_matter_parameters{iSubject});
            end
            
            % ----- Write output to (.mat) file ----- %
            filename = [save_dir '\Output_confirmatory_analysis_' seed_region(1:end-4) '_vs_' non_seed_region(1:end-4) '.mat'];
            save(filename, 'betas_seed_vs_ROI');
    end
end

end
