function white_matter_parameters = ObtainWhiteMatterPredictor(fileInfo, reslice_mask)

% Obtain white-matter parameters.
%
% Input argument(reslice_mask): Indicate whether the CSF mask to be
% eroded needs to be resliced first.
%           1 = Reslice mask
% Subfunctions: 


    % ----- Obtain white-matter indices from anatomical scans -------------
white_matter_parameters = {};
for iSubject = 1:length(fileInfo.white_matter_file_names)
    this_white_matter_file = fileInfo.white_matter_file_names{iSubject};
    
    % ----- Erode and obtain CSF anatomical file ----- %
    eroded_matrix = ErodeNuisanceVariableMask(this_white_matter_file, reslice_mask);
    
    % ----- Obtain indices of white matter voxels from segmentation ----- %
    white_matter_indices = find(eroded_matrix > 0.9);
        
    % ----- Obtain white-matter signal from functional scans ----- %
    all_functional_scans = fileInfo.functional_file_names{iSubject};
    white_matter_mean_signal = [];
    for iScan = 1:length(all_functional_scans)
        this_functional_scan = all_functional_scans{iScan};
        
        % ----- Read this functional scan ----- %
        H_Scan = spm_vol(this_functional_scan);
        D_Scan = spm_read_vols(H_Scan);
        
        % ----- Obtain mean of white matter for this subject ----- %
        white_matter_voxels = D_Scan(white_matter_indices);
        white_matter_mean = mean(white_matter_voxels);
        white_matter_mean_signal = [white_matter_mean_signal; white_matter_mean];
    end
    
    white_matter_parameters{iSubject} = white_matter_mean_signal;
end