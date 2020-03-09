function ROI_map = DefineCustomROI(input_XYZ, output_filename)

% Define costum region-of-interest.
%
% Input argument(input_XYZ): Specify the x-, y-, z-coordinates around which
% the ROI is to be defined, e.g., for a voxel in the right DLPFC:
%       input_XYZ = [1 -29 -12]
% Input argument(output_filename): Specify the name under which the output
% is to be saved, e.g.:
%       'Right_DLPFC.nii'
% Subfunctions: FileOrganizer


fileInfo = FileOrganizer(0);

% ----- Obtain the voxel size of the region-of-interest ----- %
fileID = fopen('Preprocessing_job.m', 'r');
expression_to_find = 'woptions.vox';
while ~feof(fileID)
    this_line = fgetl(fileID);
    
    % ----- Find and save the voxel size to variable ----- %
    is_it_equal = contains(this_line, expression_to_find);
    if is_it_equal == 1
        line_to_find = strsplit(this_line, ' = ');
        voxel_size = str2num(line_to_find{2});
    end
end
fclose(fileID);

% ----- Create (temporary) zeros matrix ---- %
H_Template = spm_vol(fileInfo.functional_file_names{1}{1});
[D_Template, XYZ_Template] = spm_read_vols(H_Template);
ROI_map = zeros(size(D_Template));

% ----- Specify x-, y-, z-value ranges to be included in the ROI ----- % 
x_range = [(input_XYZ(1) - voxel_size(1)) (input_XYZ(1) + voxel_size(1))];
y_range = [(input_XYZ(2) - voxel_size(2)) (input_XYZ(2) + voxel_size(2))];
z_range = [(input_XYZ(3) - voxel_size(3)) (input_XYZ(3) + voxel_size(3))];

% ----- Find indices falling within the above-defined ranges ----- %
indices_to_fill = [];
    % ----- Loop over x-coordinates within x-range ----- %
for this_x_coordinate = min(x_range):max(x_range)
        % ----- Loop over y-coordinates within y-range ----- %
    for this_y_coordinate = min(y_range):max(y_range)
            % ----- Loop over z-coordinates within z-range ----- %
        for this_z_coordinate = min(z_range):max(z_range)
            
            % ----- Collect xyz-coordinate (i.e., if it exists) ----- %
            this_xyz_index = find(XYZ_Template(1,:) == this_x_coordinate &  XYZ_Template(2,:) == this_y_coordinate & XYZ_Template(3,:) == this_z_coordinate);
            if ~isempty(this_xyz_index)
                indices_to_fill = [indices_to_fill; this_xyz_index];
            end
            
        end % End of first for-loop
    end % End of second for-loop
end % End of third for-loop

% ----- Set the collected indices to '1' in ROI map ----- %
ROI_map(indices_to_fill) = 1;

% ----- Save ROI map to file ----- %
H_Template.fname = [fileInfo.base_dir '\' output_filename];
H_Template.pinfo(1) = 1;
H_Template.dt = [16 0];
spm_write_vol(H_Template, ROI_map);
