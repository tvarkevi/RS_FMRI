function Preprocessing(base_dir, data_dir)

% Preprocessing of fMRI data
%
% Input argument(base_dir): Location of the jobfile (working directory)
% Input argument(data_dir): String name of the data directory.
% Subfunctions: inner_preprocess

% List of open inputs
%       Realign: Estimate: Session - cfg_files
%       Coregister: Estimate: Reference Image - cfg_files
%       Normalise: Estimate & Write: Image to Align - cfg_files


% ----- Selecting files ----- %
all_subjects = dir([data_dir '\xm*']);
for iSubj = 1:length(all_subjects)
    this_subj_name = all_subjects(iSubj).name;
    this_subj_path = [data_dir '\' this_subj_name];
    fprintf([num2str(iSubj) '\t' strrep(this_subj_path, '\', '\\') '\n']);
    inner_preprocess(base_dir, this_subj_path);
end

end

function inner_preprocess(base_dir, this_subj_path)

all_folders = dir([this_subj_path '\xm*']);
functional_scans = dir([this_subj_path '\' all_folders(2).name '\xm*.nii']);

inputs = cell(3,1);

% ----- Obtain functional scans ----- %
for iScan = 1:length(functional_scans)
    inputs{1}{iScan, 1} = [this_subj_path '\' all_folders(2).name '\' functional_scans(iScan).name ',1'];
end

% ----- Obtain anatomical scans ----- %
anatomical_scans = dir([this_subj_path '\' all_folders(1).name '\xm*.nii']);
inputs{2} = {[this_subj_path '\' all_folders(1).name '\' anatomical_scans(1).name ',1']};
inputs{3} = inputs{2};

jobfile = {[base_dir '\Preprocessing_job.m']};
spm_jobman('run', jobfile, inputs{:});

end
