function output_filename = LiebermanCunninghamMethod(input_filename, input_cluster_criteria)

% Conduct cluster analysis according to the Lieberman and Cunningham (2009)
% method.
%
% Input arguments (filename): Specify the name of the file for which 
% clustering is to be performed.
% Subfunctions: hiro3_islander_iterative_borg
%
% NOTE: This script requires hiro3 to be in the path in order for it to 
% work.


filename = input_filename;
cluster_criteria = input_cluster_criteria;

% ----- Read the statistical (.nii) file ----- %
H_Scan = spm_vol(filename);
D_Scan = spm_read_vols(H_Scan);

% ----- Define t cut-off value ----- %
t_cutoff = icdf('t', 1 - cluster_criteria.p_thres, cluster_criteria.df);

% ----- Perform cluster analysis ----- %
proportion_of_contact = 0.25;
[labelMap, ~] = hiro3_islander_iterative_borg(D_Scan, t_cutoff, proportion_of_contact, cluster_criteria.min_vox);

% ----- Save cluster map to (.nii) file ----- %
H_Scan.fname = ['Cluster_map_' filename];
spm_write_vol(H_Scan, labelMap);

% ----- Define output argument (filename of cluster map) ----- %
output_filename = H_Scan.fname;
