function ComputeOFCMask

% Compute probability mask for total OFC as seed region-of-interest.
% 
% Input arguments: -
% Subfunctions(ResliceTotalOFC): Reslice using SPM12


fileInfo = FileOrganizer(0);

% ----- Define base directory of probability maps ----- %
baseDir = fileInfo.base_dir;

% ----- Load headers and read data of OFC subnuclei with spm12 ----- %
OFC1Header = spm_vol([baseDir '\' 'OFC_Fo1.nii']);
OFC1Data = spm_read_vols(OFC1Header);

OFC2Header = spm_vol([baseDir '\' 'OFC_Fo2.nii']);
OFC2Data = spm_read_vols(OFC2Header);

OFC3Header = spm_vol([baseDir '\' 'OFC_Fo3.nii']);
OFC3Data = spm_read_vols(OFC3Header);

% ----- Compute and binarize total OFC seed-region ----- %
totalOFCProbabilityMap = OFC1Data + OFC2Data + OFC3Data;

% ----- Save total OFC binary mask/probability map to file ----- %
H_map = OFC1Header;
H_map.fname = [baseDir '\' 'OFC_Total.nii']; 

spm_write_vol(H_map, totalOFCProbabilityMap);

% ----- Reslice total OFC files to lower resolution ----- %
ResliceROIMask(fileInfo, H_map.fname)

end