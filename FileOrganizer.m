function fileInfo = FileOrganizer(simulate_data)

% Obtain directories and group membership for this analysis, and save into
% single struct variable
%
% Input arguments(simulate_data): Specify if real or simulated data is to
% be used, e.g.,
%       0. Use real data
%       1. Use (and create) simulated data
% Subfunctions: DataSimulation
%
% Output arguments subfields:
%   fileInfo.base_dir
%   fileInfo.data_dir
%   fileInfo.all_subjects
%   fileInfo.functional_file_names
%   fileInfo.anatomical_file_names
%   fileInfo.motion_parameters
%   fileInfo.white_matter_file_names
%   fileInfo.CSF_file_names


% ----- Define base directory ----- %
fileInfo.base_dir = 'O:\MGGZ\WO\Tim\Resting State fMRI\Analysis\Version 2 (22-01-2019)';    % Location of working directory
fileInfo.data_dir = 'G:\Resting State fMRI v.2';                                                % Location of the nifti data files

% ----- Define subject names and group assignments ----- %
filename = 'subjects.txt';

fileContents = table2cell(readtable(filename));
fileInfo.all_subjects = {};
for iColumn = 1:size(fileContents, 2)
    fileInfo.all_subjects = {fileInfo.all_subjects{:} fileContents(:, iColumn)};
end

% ----- Define directories of all functional and anatomical scans ----- %
fileInfo.functional_file_names = {};
fileInfo.anatomical_file_names = {};

    % ----- Loop over all subjects ----- %
for iSubject = 1:size(fileInfo.all_subjects{1}, 1)
    this_subject = [fileInfo.data_dir '\' fileInfo.all_subjects{1}{iSubject}];
    this_subject_group = fileInfo.all_subjects{2}{iSubject};
    
    % ----- Define subject directory ----- %
    this_subject_dir = dir([this_subject '\xm*']);
    
    % ----- Obtain all (preprocessed) fMRI scans for this subject ----- %
    this_subject_functional_dir = [this_subject_dir(2).name];
    this_subject_functional_scans = dir([this_subject '\' this_subject_functional_dir '\wxm*']);
    
    % ----- Display simulation progress; conditional statement ----- %
    if simulate_data == 1
        fprintf([num2str(iSubject) '\tSimulating data for subject: ' this_subject_dir(2).name '\t\t' datestr(now) '\n']);
    end
    
        % ----- Loop over all scans for this subject ----- %
    all_scan_names = {};
    for iScan = 1:length(this_subject_functional_scans)
        this_scan_name = [this_subject '\' this_subject_dir(2).name '\' this_subject_functional_scans(iScan).name ];
            % ----- Obtain simulated or real data ----- %
        if simulate_data == 1
            simulated_data_dir = [fileInfo.data_dir '\Data simulation'];
            if isempty(dir([simulated_data_dir '\' fileInfo.all_subjects{1}{iSubject} '\Simulated_*' num2str(iScan) '.nii']))
                this_scan_name = DataSimulation(iScan, this_scan_name, this_subject_group, fileInfo.base_dir, simulated_data_dir);
            else
                this_subject_simulated_dir = [simulated_data_dir '\' fileInfo.all_subjects{1}{iSubject}];
                this_subject_simulated_scans = dir([this_subject_simulated_dir '\*xm*.nii' ]);
                this_scan_name = [simulated_data_dir '\' fileInfo.all_subjects{1}{iSubject} '\' this_subject_simulated_scans(iScan).name];
            end
            all_scan_names{iScan} = this_scan_name;
        elseif simulate_data == 0
            all_scan_names{iScan} = this_scan_name;
        end
    end
    fileInfo.functional_file_names{iSubject} = all_scan_names;
    
    % ----- Obtain motion parameters for this subject ----- %
    this_subject_motion_parameters_dir = dir([this_subject '\' this_subject_functional_dir '\rp*.txt']);
    fileID = fopen([this_subject '\' this_subject_functional_dir '\' this_subject_motion_parameters_dir.name], 'r');    
    all_parameters_cell = textscan(fileID, '%s');
    all_parameters = [];
    for iParameter = 1:length(all_parameters_cell{1})
        all_parameters = [all_parameters; str2double(all_parameters_cell{1}{iParameter})];
    end
    fileInfo.motion_parameters{iSubject} = [all_parameters(1:6:end) ...
        all_parameters(2:6:end) all_parameters(3:6:end) ...
        all_parameters(4:6:end) all_parameters(5:6:end) ...
        all_parameters(6:6:end)];
    fclose(fileID);
    
    % ----- Obtain (preprocessed) anatomical scan for this subject ----- %
    this_subject_anatomical_dir = [this_subject_dir(1).name];
    this_subject_anatomical_scan = dir([this_subject '\' this_subject_anatomical_dir '\wxm*.nii']);
    fileInfo.anatomical_file_names{iSubject} = [this_subject '\' this_subject_anatomical_dir '\' this_subject_anatomical_scan.name];
       
    % ----- Obtain white matter parameters -----%
    this_subject_white_matter_dir = [this_subject '\' this_subject_anatomical_dir]; 
    this_subject_white_matter_scan = dir([this_subject_white_matter_dir '\c2*']); 
    fileInfo.white_matter_file_names{iSubject} = [this_subject_white_matter_dir '\' this_subject_white_matter_scan.name];
    
    % ----- Obtain CSF parameters -----%
    this_subject_csf_dir = [this_subject '\' this_subject_anatomical_dir];
    this_subject_csf_scan = dir([this_subject_csf_dir '\c3*']);
    fileInfo.CSF_file_names{iSubject} = [this_subject_csf_dir '\' this_subject_csf_scan.name];
 end
