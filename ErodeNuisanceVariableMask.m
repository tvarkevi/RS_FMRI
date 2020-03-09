function eroded_matrix = ErodeNuisanceVariableMask(input_filename, reslice_mask)

% Erosion of segmented white matter or CSF anatomical data.
%
% Input argument (input_filename): Specify the name of the file to be
% eroded, e.g.:
%       'H:\Psychopatie project\Nifti databestanden\xm13101101\xm13101101_4_1\c2wxm13101101_4_1-0001.nii'
% Input argument (reslice_mask): Indicate whether the nuisance variable map
% needs to be resliced, prior to erosion:
%       1. Reslice the nuisance variable map prior to erosion
% Subfunctions: FileOrganizer, ResliceNuisanceVariableMask


fileInfo = FileOrganizer(0);

% ----- Reslice nuisance variable mask ----- %
if reslice_mask == 1
    ResliceNuisanceVariableMask(fileInfo, input_filename)
end

% ----- Rename input filename ----- %
filenameParts = strsplit(input_filename, '\');
filename = [];
for iParts = 1:length(filenameParts) - 1
    filename = [filename filenameParts{iParts} '\'];
end
filename = [filename 'r' filenameParts{end}];

% ----- Obtain (meta)data of the nuisance mask itself ----- %
H_Scan = spm_vol(filename);
D_Scan = spm_read_vols(H_Scan);
D_Scan(D_Scan > 0.8) = 1;
D_Scan(D_Scan <= 0.8) = 0;

% ----- Create empty matrices with size identical to original data ----- %
horizontal_matrix = zeros(size(D_Scan));
sagittal_matrix = zeros(size(D_Scan));
coronal_matrix = zeros(size(D_Scan));

[nRows, nCols, nSlices] = size(D_Scan);
    % ----- Erode in horizontal plane ----- %
for iHorizontalSlice = 2:(nSlices-1)
    this_horizontal_slice = D_Scan(:,:,iHorizontalSlice);
    
    % ----- Obtain slices above and below the current slice ----- %
    slice_below = D_Scan(:,:,(iHorizontalSlice - 1));
    slice_above = D_Scan(:,:,(iHorizontalSlice + 1));
    
    % ----- Multiply with slices above and below current slice ----- %
    this_horizontal_slice = this_horizontal_slice .* slice_below .* slice_above;
    
    % ----- Replace current slice with resampled slice ----- %
    horizontal_matrix(:,:,iHorizontalSlice) = this_horizontal_slice;
end

    % ----- Erode in coronal plane ----- %
for iCoronalSlice = 2:(nCols-1)
    this_coronal_slice = D_Scan(:,iCoronalSlice,:);
    
    % ----- Reshape current slice ----- %
    this_coronal_slice = reshape(this_coronal_slice,[nRows,nSlices]);
    
    % ----- Obtain (and reshape) slices behind the current slice ----- %
    slice_behind = D_Scan(:,(iCoronalSlice - 1),:);
    slice_behind = reshape(slice_behind,[nRows,nSlices]);
    
    % ----- Obtain (and reshape) slices in front the current slice ----- %
    slice_front = D_Scan(:,(iCoronalSlice + 1),:);
    slice_front = reshape(slice_front,[nRows,nSlices]);
    
    % ----- Multiply current slice with slice in front and behind ----- %
    this_coronal_slice = this_coronal_slice .* slice_behind .* slice_front;
    
    % ----- Replace current slice with resampled slice ----- %
    coronal_matrix(:,iCoronalSlice,:) = this_coronal_slice;
end

    % ----- Erode in sagittal plane ----- %
for iSagittal = 2:(nRows-1)
    this_sagittal_slice = D_Scan(iSagittal,:,:);

    % ----- Reshape current slice ----- %
    this_sagittal_slice = reshape(this_sagittal_slice,[nCols,nSlices]);
    
    % ----- Obtain (and reshape) slices left from current slice ----- %
    slice_left = D_Scan((iSagittal - 1),:,:);
    slice_left = reshape(slice_left,[nCols,nSlices]);
    
    % ----- Obtain (and reshape) slices right from current slice ----- %
    slice_right = D_Scan((iSagittal + 1),:,:);
    slice_right = reshape(slice_right,[nCols,nSlices]);
    
    % ----- Multiply current slice with slices left and right ----- %
    this_sagittal_slice = this_sagittal_slice .* slice_left .* slice_right;
    
    % ----- Replace current slice with resampled slice ----- %
    sagittal_matrix(iSagittal,:,:) = this_sagittal_slice;
end

% ----- Multiply horizontal, coronal, and sagittal matrices ----- %
eroded_matrix = coronal_matrix .* sagittal_matrix .* horizontal_matrix;

% ----- Save eroded mask to (.nii) file ----- %
filenameParts = strsplit(filename, '\');
H_Scan.fname = [];
for iParts = 1:length(filenameParts) - 1
    H_Scan.fname = [H_Scan.fname filenameParts{iParts} '\'];
end
H_Scan.fname = [H_Scan.fname 'Eroded_' filenameParts{end}];
% H.pinfo(1) = 1;
spm_write_vol(H_Scan, eroded_matrix);

end

function ResliceNuisanceVariableMask(fileInfo, input_filename)

inputs{1} = {fileInfo.functional_file_names{1}{1}};
inputs{2} = {[input_filename ',1']};

jobfile = {[fileInfo.base_dir '\ResliceROIMask_job.m']};
spm_jobman('run', jobfile, inputs{:});

end