function ComputeAmygdalaMask

% Compute probability mask for total amygdala as seed region-of-interest.
% 
% Input arguments: -
% Subfunctions(ResliceROIMask): Reslice using SPM12


fileInfo = FileOrganizer(0);

% ----- Define base directory of probability maps ----- %
baseDir = fileInfo.base_dir;

% ----- Load headers and read data of amygdala subnuclei with spm12 ----- %
basolateralAmygdalaHeader = spm_vol([baseDir '\' 'Amygdala_LB.nii']);
basolateralAmygdalaData = spm_read_vols(basolateralAmygdalaHeader);

centromedialAmygdalaHeader = spm_vol([baseDir '\' 'Amygdala_CM.nii']);
centromedialAmygdalaData = spm_read_vols(centromedialAmygdalaHeader);

superficialAmygdalaHeader = spm_vol([baseDir '\' 'Amygdala_SF.nii']);
superficialAmygdalaData = spm_read_vols(superficialAmygdalaHeader);

% ----- Compute and binarize total amygdala seed-region ----- %
totalAmygdalaProbabilityMap = basolateralAmygdalaData + centromedialAmygdalaData + superficialAmygdalaData;

% ----- Save total amygdala binary mask/probability map to file ----- %
H_map = basolateralAmygdalaHeader;
H_map.fname = [baseDir '\' 'Amygdala_Total.nii']; 

spm_write_vol(H_map, totalAmygdalaProbabilityMap);

% ----- Reslice total amygdala files to lower resolution ----- %
ResliceROIMask(fileInfo, H_map.fname)

end