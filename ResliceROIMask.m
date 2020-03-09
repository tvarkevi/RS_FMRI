function ResliceROIMask(fileInfo, input_filename)

% Reslice the ROI mask.
%
% Input argument(fileInfo): The functional filenames are used.
% Subfunctions: ResliceROIMask_job


inputs{1} = {fileInfo.functional_file_names{1}{1}};
inputs{2} = {[input_filename ',1']};

jobfile = {[fileInfo.base_dir '\ResliceROIMask_job.m']};
spm_jobman('run', jobfile, inputs{:});

end