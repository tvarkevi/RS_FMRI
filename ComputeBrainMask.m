function ComputeBrainMask(fileInfo, reslice_grey_matter_mask)

% Compute brain mask for each subject
%
% Input argument: fileInfo
% Input argument(reslice_grey_matter_mask): Reslice the segmented grey
% matter files, e.g.:
%       1: Reslice the grey matter masks
% Subfunctions: -


% ----- Obtain file directories of segmentations ----- %
for iSubject = 1:length(fileInfo.anatomical_file_names)
    this_subject_anatomical_dir = fileparts(fileInfo.anatomical_file_names{iSubject});
    all_anatomical_files = dir([this_subject_anatomical_dir '\rc*.nii']);
    
    % ----- Display progress (this subject) ----- %
    disp([num2str(iSubject) 9 this_subject_anatomical_dir]); % The ASCI character for tab is 9.
    
    % ----- Reslice grey matter mask ----- %
    if reslice_grey_matter_mask == 1
        file_parts = strsplit(fileInfo.anatomical_file_names{iSubject}, '\');
        grey_matter_file = dir([strjoin(file_parts(1:end-1), '\') '\c1*.nii']);
        ResliceGreyMatterMask(fileInfo, [strjoin(file_parts(1:end-1), '\') '\' grey_matter_file.name])
    end
    
    % ----- Store segmentation directories to variable ----- %
    grey_matter_files{iSubject} = [this_subject_anatomical_dir '\' all_anatomical_files(1).name];
    white_matter_files{iSubject} = [this_subject_anatomical_dir '\' all_anatomical_files(2).name];
    CSF_files{iSubject} = [this_subject_anatomical_dir '\' all_anatomical_files(3).name];
    
    % ----- Read subject-level segmentation masks ----- %
    H_Brain = spm_vol(grey_matter_files{iSubject});
    D_Brain = spm_read_vols(H_Brain);
    H_White = spm_vol(white_matter_files{iSubject});
    D_White = spm_read_vols(H_White);
    H_CSF = spm_vol(CSF_files{iSubject});
    D_CSF = spm_read_vols(H_CSF);
    
    % ----- Compute brain mask for this subject ----- %
    brain_mask_raw = D_Brain + D_White + D_CSF;
    find_brain_voxels = find(brain_mask_raw > 0);
    brain_mask = zeros(size(brain_mask_raw));
    brain_mask(find_brain_voxels) = 1;
    
    [nRows, nCols, nSlices] = size(brain_mask);
    for iVoxel = 2:((nRows * nCols * nSlices) - 1)
        previous_voxel = brain_mask(iVoxel - 1);
        next_voxel = brain_mask(iVoxel + 1);
        % ----- Set voxel to one if adjacent voxels are also one ----- %
        if previous_voxel == 1 && next_voxel == 1
            brain_mask(iVoxel) = 1;
        else
            continue
        end
    end
    
    % ----- Save brain mask to (.nii) file ----- %
    H_Brain.fname = [this_subject_anatomical_dir '\Brain_Mask.nii'];
    spm_write_vol(H_Brain, brain_mask);
end

end

function ResliceGreyMatterMask(fileInfo, input_filename)

inputs{1} = {fileInfo.functional_file_names{1}{1}};
inputs{2} = {[input_filename ',1']};

jobfile = {[fileInfo.base_dir '\ResliceROIMask_job.m']};
spm_jobman('run', jobfile, inputs{:});

end
