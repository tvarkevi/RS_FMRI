function ComputeAICMask

% Compute probability mask for AIC as seed region-of-interest.
% 
% Input arguments: -
% Subfunctions(ResliceTotalOFC): Reslice using SPM12


fileInfo = FileOrganizer(0);

% ----- Define base directory of probability maps ----- %
baseDir = fileInfo.base_dir;

% ----- Reslice the AICHA atlas to the native resolution ----- %
ResliceROIMask(fileInfo, [baseDir '\' 'AICHAmc.nii']);

% ----- Load headers and read data of the resliced AICHA atlas ----- %
H_rAICHA = spm_vol([baseDir '\' 'rAICHAmc.nii']);
D_rAICHA = spm_read_vols(H_rAICHA);

% ----- Find anterior insular indices ----- %
allAICIndices = find(D_rAICHA == 74 | D_rAICHA == 75 | D_rAICHA == 76 | D_rAICHA == 77 | D_rAICHA == 78);

% ----- Compute AIC mask ----- %
D_Map = zeros(size(D_rAICHA));
D_Map(allAICIndices) = 1;

% ----- Save new mask to file ----- %
H_map = H_rAICHA;
H_map.fname = [baseDir '\' 'rAIC.nii'];

spm_write_vol(H_map, D_Map);
