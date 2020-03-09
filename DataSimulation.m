function this_scan_name = DataSimulation(change_factor, input_filename, group_membership, base_dir, simulated_data_dir)

% Simulate data for software testing.
%
% Input argument(change_factor): Specify by which factor increase/decrease
% amount will be multiplied, and added to the actual voxel value.
% Input argument(input_filename): Specify the input filename, e.g.,
%       'H:\Psychopatie project\Nifti databestanden\xm13101101\xm13101101_5_1\wxm13101101_5_1-0001.nii'
% Input argument(group_membership): Indicate to which group the subject
% belongs, e.g.,
%       0. Low psychopathy/Control subject
%       1. High psychopathy/Aggression subject
% Input argument(base_dir): Specify working directory; i.e., location of 
% the (seed) region-of-interest maps. E.g.,
%       'O:\MGGZ\WO\Nadia\Psychopatieproject - Scripts'
% Input argument(data_dir): Specify directory where data is (to be) stored,
% e.g.,
%       'H:\Psychopatie project\Nifti databestanden'
% Subfunction: -


% ----- Specify factor to base simulated correlation on ----- %
baseline_voxel_value = 1;
absolute_change = 8 / 320;

% ----- Read scan on which simulation is to be based ----- %
H_Scan = spm_vol(input_filename);
D_Scan = spm_read_vols(H_Scan);
simulated_data = randn(size(D_Scan));

% ----- Find indices of voxels outside brain ----- %
file_parts = strsplit(input_filename, '\');
this_subject_dir = strjoin(file_parts(1:end-2), '\');
find_anatomical_dir = dir([this_subject_dir '\xm*']);
H_Mask = spm_vol([this_subject_dir '\' find_anatomical_dir(1).name '\Brain_Mask.nii']);
D_Mask = spm_read_vols(H_Mask);
non_brain_indices = find(D_Mask < 1);

% ----- Set non-brain voxels to NaN ----- %
simulated_data(non_brain_indices) = NaN;

% ----- Find indices of seed-region ----- %
seed_region = 'rAmygdala_Total.nii';
H_Seed = spm_vol([base_dir '\' seed_region]);
D_Seed = spm_read_vols(H_Seed);
seed_indices = find(D_Seed > 0.2);

% ----- Apply increase factor to seed-region (amygdala) voxels ----- %
% simulated_data(seed_indices) = simulated_data(seed_indices) + (change_factor * absolute_change);
simulated_data(seed_indices) = baseline_voxel_value + (change_factor * absolute_change);

% ----- Find indices of (non-seed) ROI ----- %
non_seed_ROI = 'rOFC_Total.nii';
H_ROI = spm_vol([base_dir '\' non_seed_ROI]);
D_ROI = spm_read_vols(H_ROI);
non_seed_ROI_indices = find(D_ROI > 0);

% ----- Apply decrease factor to (non-seed) ROI voxels ---- %
if group_membership == 0
    simulated_data(non_seed_ROI_indices) = baseline_voxel_value - (rand / 5) - (change_factor * absolute_change);
elseif group_membership == 1
    simulated_data(non_seed_ROI_indices) = baseline_voxel_value - (change_factor * absolute_change);
end

% ----- Create data simulation (sub)folder(s) if not existent ----- %
if ~exist(simulated_data_dir)
    mkdir(simulated_data_dir);
end
new_subject_directory = [simulated_data_dir '\' file_parts{end-2}];
if ~exist(new_subject_directory)
    mkdir(new_subject_directory);
end

% ----- Write simulated data/scan to (.nii) file ----- %
H_Scan.fname = [new_subject_directory '\Simulated_' file_parts{end}];
this_scan_name = H_Scan.fname;
H_Scan.pinfo(1) = 1;
H_Scan.dt = [16 0];
spm_write_vol(H_Scan, simulated_data);
