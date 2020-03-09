function CSF_parameters = ObtainCSFPredictor(fileInfo, reslice_mask)

% Obtain cerebrospinal fluid parameters.
%
% Input argument(reslice_mask): Indicate whether the CSF mask to be
% eroded needs to be resliced first.
%           1 = Reslice mask
% Subfunctions: -

    % ----- Obtain CSF indices from anatomical scans ----------------------
CSF_parameters = {};
for iSubject = 1:length(fileInfo.CSF_file_names)
    this_CSF_file = fileInfo.CSF_file_names{iSubject};
    
    % ----- Erode and obtain CSF anatomical file ----- %
    eroded_matrix = ErodeNuisanceVariableMask(this_CSF_file, reslice_mask);
    
    % ----- Obtain indices of CSF voxels from segmentation ----- %
    CSF_indices = find(eroded_matrix > 0.9);
        
    % ----- Obtain CSF signal from functional scans -----------------------
    all_functional_scans = fileInfo.functional_file_names{iSubject};
    CSF_mean_signal = [];
    for iScan = 1:length(all_functional_scans)
        this_functional_scan = all_functional_scans{iScan};
        
        % ----- Read this functional scan ----- %
        H_Scan = spm_vol(this_functional_scan);
        D_Scan = spm_read_vols(H_Scan);
        
        % ----- Obtain mean of CSF for this subject ----- %
        CSF_voxels = D_Scan(CSF_indices);
        CSF_mean = mean(CSF_voxels);
        CSF_mean_signal = [CSF_mean_signal; CSF_mean];
    end
    
    CSF_parameters{iSubject} = CSF_mean_signal;
end