function Segmentation(base_dir, data_dir)

% Normalize and Segment the anatomical scans into grey-matter (c1),
% white-matter (c2), and CSF (c3).
%
% Input arguments: -
% Subfunctions: Segmentation_job


% ----- Obtain all subject directories ----- %
all_subjects = dir([data_dir '\xm*']);

    % ----- Loop over all subjects in all_subjects ----- %
for iSubject = 1:length(all_subjects)
    this_subject_dir = [data_dir '\' all_subjects(iSubject).name];
        
    % ----- Obtain all subfolders for this subject ----- %
    this_subject_all_folders = dir([this_subject_dir '\xm*']);
    fprintf([num2str(iSubject) '\t' strrep(this_subject_dir, '\', '\\') '\n']);
    
    % ----- Obtain T1 anatomical scan for this subject ----- %
    this_subject_anatomy_dir = [this_subject_dir '\' this_subject_all_folders(1).name];
    this_subject_anatomy_scan = dir([this_subject_anatomy_dir '\xm*.nii']);    
    
    % ----- Specify inputs for spm segmentation for this subject ----- %
    inputs{1} = {[this_subject_anatomy_dir '\' this_subject_anatomy_scan.name ',1']};
    inputs{2} = inputs{1};
    
    % ----- Conduct spm segmentation for this subject ----- %
    jobfile = {[base_dir '\Segmentation_job.m']}; 
    spm_jobman('run', jobfile, inputs{:});
end
