# RS_FMRI
Seed-based functional connectivity analysis of resting-state fMRI

Preliminary note: The resting-state seed-based functional connectivity procedures described in this document were developed on a Windows 7 workstation in Matlab R2016b. The following auxiliary software and toolboxes are required to be installed and added to the Matlab path in order for this analysis suite to work: SPM12, SPM Anatomy toolbox, r2agui, hiro3.

## 1. Data preprocessing and preparation

The analyses described in this document assume you have the raw (PAR/REC) MR data files converted into NIFTI format (e.g., by using the r2agui toolbox; see http://r2agui.sourceforge.net/) and stored (list-wise) in a directory called data_dir.

### 1.1 Data management

Many of the functions described in this manual call a script called FileOrganizer.m. This function contains all hard-coded filenames and full paths that the main analysis pipeline uses. Hence, in order to make use of the resting-state seed-based connectivity software described here, it is first necessary to change lines 25 (i.e., working/base directory), 26 (i.e., data directory), 29 (i.e., the name of the text file that contains the list of subject names), 47, 51, 66, and 92 (i.e., naming conventions of original or resliced scan files) in FileOrganizer.m so that they match the file infrastructure and naming conventions of your specific system.

Note: The software assumes the scripts, subjects.txt, and ROI maps are stored in the directory specified by base_dir, with the data stored separately (and list-wise) in the directory specified by data_dir.

To run the FileOrganizer.m script separately from the main analysis suite, enter the following lines of code in the Matlab command window:

```
simulate_data = 0
fileInfo = FileOrganizer(simulate_data);
```

From this point onwards, it is assumed that fileInfo is a variable stored in workspace, with no alterations being performed on it, unless specifically specified in this manual.

### 1.2 Preprocessing and segmentation

As a first step, the data need to be preprocessed and (the anatomical data) segmented into grey matter, white matter, and cerebrospinal fluid maps. This is done using the Preprocessing.m and Segmentation.m scripts.
For the present purposes, the raw (functional) data are subjected to the following preprocessing steps:
1. Motion correction (realignment)
2. Co-registration to the anatomical image
3. Normalization to standard MNI space
To conduct these 3 preprocessing steps, enter the following lines of code in the Matlab command window:

```
Preprocessing(fileInfo.base_dir, fileInfo.data_dir);
```

Note: In order to successfully run the Preprocessing.m script, it is necessary to first change line 25 of the preprocessing jobfile, Preprocessing_job.m, into the full path of the TPM.nii file located on your device (usually in the spm12/tpm subdirectory).

To segment the anatomical data into grey matter, white matter, and cerebrospinal fluid maps, enter the following line of code into the command window:

```
Segmentation(fileInfo.base_dir, fileInfo.data_dir);
```

Note: In order to successfully run the Segmentation.m script, it is first necessary to change lines 10, 24, 28, 32, 36, 40, and 44 of the segmentation jobfile, Segmentation_job.m, into the full path of the TPM.nii file located on your device (usually in the spm12/tpm subdirectory).

### 1.3 Brain mask

As a last preparatory step, compute a mask that contains only voxels within the brain, for each of the participants specified in the subjects.txt file. This is done via the ComputeBrainMask.m script. Enter the following lines of code in the Matlab command window to execute the function:

```
reslice_grey_matter_mask = 1;
ComputeBrainMask(fileInfo, reslice_grey_matter_mask);
```

## 2. Defining the regions-of-interest

2.1 Amygdala (seed) regions

Reslice the basolateral and centromedial amygdala (probability, or p) maps (respectively). The basolateral amygdala p-map is resliced (via the ResliceROIMask.m script) by entering the following lines of code in the Matlab command window:

input_filename = 'Amygdala_LB.nii';
ResliceROIMask(fileInfo, input_filename);
	
Next, reslice the centromedial amygdala p-map by entering the following lines of code in the command window (i.e., assuming fileInfo is already defined; see above):

input_filename = 'Amygdala_CM.nii';
ResliceROIMask(fileInfo, input_filename);

Note: It is also possible to create a p-map of the total amygdala as seed-region, This can be done using the ComputeAmygdalaMask.m file by entering the following line of code in the command window:

ComputeAmygdalaMask;

This function sums the ‘Amygdala_BL.nii’, ‘Amygdala_CM.nii’, and ‘Amygdala_SF.nii’ p-maps into the conjoined output files ‘Amygala_Total.nii’ and ‘Amygala_Total.nii’, representing the unprocessed and resliced total amygdala p-maps, respectively.

2.2 Anterior cingulate cortex (ACC)

Reslice the ACC p-map by entering the following lines of code into the command window:

input_filename = 'Cingulum_33.nii';
ResliceROIMask(fileInfo, input_filename);

2.3 Anterior insular cortex (AIC)
The Anatomy toolbox of spm12 does not contain a p-map of the AIC. Hence, the AIC is defined as seed-region by using the AICHA connectivity atlas of Joliot et al. (2015). To compute a (resliced) binary map of the AIC as region-of-interest, enter the following line of code in the command window:

ComputeAICMask;

2.4 Orbitofrontal cortex (OFC)

An ROI-wise (group) analysis cannot be conducted without first computing the (non-seed) ROIs themselves. To compute a (resliced) p-map of the total OFC as region-of-interest, enter the following line of code in the command window:

ComputeOFCMask;

This ComputeOFCMask.m script generates a (resliced) p-map of the OFC as the sum of the ‘OFC_Fo1.nii’, ‘OFC_Fo2.nii’, and ‘OFC_Fo3.nii’  data maps from the Anatomy toolbox of spm12. The output is stored into the ‘OFC_Total.nii’ and ‘rOFC_Total.nii’ files, which represent the unprocessed and resliced total OFC p-maps, respectively.

2.5 Periaqueductal grey (PAG)

The Anatomy toolbox of spm12 does not contain a p-map of the PAG. Therefore, to use the PAG as region-of-interest, a set of coordinates derived from a recent meta-analysis by Linnman et al. (2012) are used: x = 1, y = -29, z = -12. A binary spherical ROI is built around these central coordinates via the DefineCustomROI.m script. Enter the following lines of code in the command window:

input_XYZ = [1 -29 -12];
output_filename = 'PAG.nii';
ROI_map = DefineCustomROI(input_XYZ, output_filename);

3. Confirmatory ROI analysis

The confirmatory ROI analyses are conducted in a 3-step procedure:
1. The time-course of the (seed-)ROIs and nuisance variables are extracted from the preprocessed functional data (see section 3.2.3).
2. Subject-level general linear models are fitted with the (filtered) time-course of the (non-seed) ROI as outcome variable, the (filtered) time-course of the seed-region as predictor-of-interest, and the six motion parameters, along with the (filtered) time-courses of the white matter, CSF, and global mean signal as nuisance variables (see section 3.2). An optional additional brain region can be entered as extra nuisance variable.
3. Independent samples t-tests are used to determine whether group differences in (mean) functional connectivity are present (see section 3.3).
These steps (as well as several auxiliary procedures) are performed via the GroupLevelAnalysis.m and SubjectLevelAnalysis.m scripts.

3.1 White matter and cerebrospinal fluid (CSF) (Optional)

The mean signal of the white matter and CSF are used as nuisance variables in the subject-level GLMs. The software pipeline can compute and utilize these nuisance variables automatically, without user-intervention via the ObtainWhiteMatterPredictor.m and ObtainCSFPredictor.m scripts, respectively.

Alternatively (but optionally), the white matter and CSF time-courses can be extracted from the data and stored as .mat files by putting lines 25 and 26 of SubjectLevelAnalysis.m in comment, than de-commenting lines 27 and 28. Next, enter the following lines of code in the command window:

reslice_mask = 1;
white_matter_parameters = ObtainWhiteMatterPredictor(fileInfo, reslice_mask);
save('white_matter_parameters.mat', 'white_matter_parameters');
CSF_parameters = ObtainCSFPredictor(fileInfo, reslice_mask);
save('CSF_parameters.mat', 'CSF_parameters');

By using this code, the white matter and CSF nuisance parameters will be stored in the white_matter_parameters.mat and CSF_parameters.mat Matlab files, respectively.

Note: the mean white matter and CSF signals are computed from the corresponding segmented anatomical data files, eroded by one voxel in each axis via the ErodeNuisanceVariableMask.m script.

3.2 Subject-level statistical analyses

The subject-level statistical analysis is conducted in a multi-step procedure:
- The SubjectLevelAnalysis.m script loops over all subjects and stores the GLM statistics in a .mat file in the functional (subjects) data directories (see section 3.2.1).
- The ROIWiseGLM.m script conducts the actual seed-to-ROI-wise GLM analysis for the current subject (see section 3.2.2).
- The ObtainROIPredictor.m script extracts the time-course of the seed- and non-seed ROIs (see section 3.2.3.

3.2.1 Subject-level analyses

The subject-level confirmatory analyses are governed by a script called SubjectLevelAnalysis.m. This script batches the seed-to-ROI-wise GLM analyses performed by the ROIWiseGLM.m script. It then saves the resulting beta-coefficients of the predictor of interest (i.e., the seed-region) as a .mat file in the subject directories. 

Optional: It is possible to enter a second brain region as nuisance variable in the seed-to-ROI-wise GLM analyses, by entering an optional additional input argument for the SubjectLevelAnalysis.m and ROIWiseGLM.m functions (i.e., via the definition statement varargin in Matlab). This approach is one of two methods to control for partial volume sharing with a region-of-no-interest that is provided by the analysis pipeline (see section 3.2.3 for the second method; see also Varkevisser, Gladwin, Heesink, van Honk, & Geuze, 2017). 

Optional: As already alluded to in section 3.1, the subject-level GLM analyses are conducted with the (filtered) white matter and CSF parameters as nuisance variables. The SubjectLevelAnalysis.m script can either extract these nuisance variables automatically (if lines 26 and 27 are in comment, but lines 24 and 25 are de-commented) or load these parameters from pre-computed .mat files (if lines 26 and 27 are de-commented, but lines 24 and 25 are in comment). If the former of these options is selected, then the ObtainWhiteMatterPredictor.m and ObtainCSFPredictor.m scripts are executed automatically by the SubjectLevelAnalysis.m function (see also section 3.1).

3.2.1.1 Basolateral amygdala vs. OFC

To execute the subject-level confirmatory GLM procedure with the basolateral amygdala as seed-region, the centromedial amygdala as optional nuisance variable, and the OFC as non-seed ROI, enter the following lines of code in the command window:

seed_region = 'rAmygdala_LB.nii';
nuisance_region = 'rAmygdala_CM.nii';
non_seed_region = 'rOFC_Total.nii';
analysis_type = 2;
SubjectLevelAnalysis(seed_region, non_seed_region, analysis_type, fileInfo, nuisance_region);

3.2.1.2 Centromedial amygdala vs. PAG

To execute the subject-level confirmatory GLM procedure with the centromedial amygdala as seed-region, the basolateral amygdala as optional nuisance variable, and the PAG as non-seed ROI, enter the following lines of code in the command window:

seed_region = 'rAmygdala_CM.nii';
nuisance_region = 'rAmygdala_LB.nii';
non_seed_region = 'PAG.nii';
analysis_type = 2;
SubjectLevelAnalysis(seed_region, non_seed_region, analysis_type, fileInfo, nuisance_region);

3.2.2 Seed-to-ROI-wise GLM analysis

The seed-to-ROI-wise GLM analysis is conducted using the ROIWiseGLM.m script, which requires a subject number (corresponding to the line of the current subject in the subjects text file) as input argument. This function first extracts the raw time-course of the (mean) seed-region signal via the ObtainROIPredictor.m script. It then filters the signal of the seed-region predictor of interest, as well as that of the nuisance variables, by using the ezfilt.m sub-function. It also computes (and filters) the global mean signal as nuisance variable. The time-course of the non-seed ROI is then extracted using the ObtainROIPredictor.m script, and again filtered with the ezfilt.m script. These filtered signals are then used to conduct the following GLM:

Non-Seed ROI = Intercept + Seed + Global Mean + Motion Parameters + White Matter + CSF + error

This GLM is computed for all (hemispheric) combinations of seed/non-seed ROIs.

To execute the ROIWiseGLM.m function for a single subject by hand, enter the following code in the command window (iSubject = 32 is provided as an example):

seed_region = 'rAmygdala_Total.nii';
non_seed_region = 'rOFC_Total.nii';
iSubject = 32;
load('CSF_parameters.mat');
CSF_parameters = CSF_parameters{iSubject};
load('white_matter_parameters.mat');
white_matter_parameters = white_matter_parameters{iSubject};
betas_seed_vs_ROI = ROIWiseGLM(seed_region, non_seed_region, iSubject, fileInfo, CSF_parameters, white_matter_parameters);

Note: The ezfilt.m script is a function created by Thomas Gladwin (https://www.tegladwin.com). It high- (output argument =  fast) and low-pass (output argument = slow) filters the variable 'signal' at the specified critical value (cutFreq), using a frequency (Fs) defined as 1/TR (i.e., the repetition time of the scans). For the present purposes, the signal will be band-pass filtered at 0.01-0.08 Hz. An example of how ezfilt.m can be used for a given signal is provided below:

TR = 1.6;
[~, signal] = ezfilt(signal, 1/TR, 0.01);
[signal, ~] = ezfilt(signal, 1/TR, 0.08);

3.2.3 Extraction of the predictors-of-interest

As mentioned above, the raw time-series of the predictor variables of interest are extracted via the ObtainROIPredictor.m script. This function reads the probability/binary atlas maps of the (seed-) regions-of-interest and extracts the BOLD signal of this ROI for each hemisphere separately. 

Optional: The ObtainROIPredictor.m script provides a second option to control for the effects of partial volume sharing via the definition statement varargin (see also section 3.2.1). More specifically, by entering an optional additional nuisance region as input argument, the function is able to control for partial volume sharing by selecting only those voxels for which the rule seed_probability_map > nuisance_probability_map applies (see Varkevisser, Gladwin, Heesink, van Honk, & Geuze, 2017). 

To execute the ObtainROIPredictor.m function for a single subject by hand, enter the following code in the command window (iSubject = 19 is provided as an example):

iSubject = 19;
nScans = length(fileInfo.functional_file_names{iSubject});
for iScan = 1:nScans
    H_Scan = spm_vol(fileInfo.functional_file_names{iSubject}{iScan});
    D_Scan = spm_read_vols(H_Scan);
    all_subject_scans{iScan} = D_Scan;
end	
region_of_interest = 'rAmygdala_LB.nii';
nuisance_region = 'rAmygdala_CM.nii';
[predictor_vector_right_ROI, predictor_vector_left_ROI] = ObtainROIPredictor(fileInfo, region_of_interest, all_subject_scans, nuisance_region);

3.3 Group-level statistical analyses

The group-level statistical analyses are governed by a script called GroupLevelAnalysis.m. In this function, the confirmatory group ROI analyses are conducted for all left/right combinations of seed-region vs. non-seed ROI, via independent samples t-tests. The statistical output data are saved in a comma-separated file (CSV) in the format: ‘Results_confirmatory_group_analyses’. 

3.3.1 Basolateral amygdala vs. OFC

To execute a confirmatory group ROI analysis with the basolateral amygdala as seed-region, and the OFC as non-seed ROI, enter the following lines of code in the command window:

seed_region = 'rAmygdala_LB.nii';
non_seed_region = 'rOFC_Total.nii';
analysis_type = 2;
GroupLevelAnalysis(fileInfo, seed_region, non_seed_region, analysis_type);

3.3.1 Centromedial amygdala vs. PAG

To execute a confirmatory group ROI analysis with the centromedial amygdala as seed-region, and the PAG as non-seed ROI, enter the following lines of code in the command window:

seed_region = 'rAmygdala_CM.nii';
non_seed_region = 'PAG.nii';
analysis_type = 2;
GroupLevelAnalysis(fileInfo, seed_region, non_seed_region, analysis_type);

4. Exploratory whole-brain analysis

The exploratory whole-brain analyses are conducted in a 4-step procedure:
1. The time-course of the (seed-)ROIs and nuisance variables are extracted from the preprocessed functional data (see section 4.2.3).
2. Subject-level beta-weight maps are created via voxel-wise GLM analyses that use the (filtered) time-course of each voxel as outcome variable, the (filtered) time-course of the seed-region as predictor-of-interest, and the six motion parameters, along with the (filtered) time-courses of the white matter, CSF, and global mean signal as nuisance variables (see section 4.2). An optional additional brain region can also be entered as extra nuisance variable.
3. Group-level statistical t-maps are created both for within and between subjects (see section 4.3). 
4. Cluster analyses are conducted to determine whether clustering of voxel-wise significance occurs in the between-groups t-map (see section 4.4)
These steps (as well as several auxiliary procedures) are performed by the GroupLevelAnalysis.m and SubjectLevelAnalysis.m scripts.

4.1 White matter and cerebrospinal fluid (CSF) (Optional)

The mean signal of the white matter and CSF are used as nuisance variables in the subject-level GLMs. The software pipeline can compute and utilize these nuisance variables automatically, without user-intervention via the ObtainWhiteMatterPredictor.m and ObtainCSFPredictor.m scripts, respectively.

Alternatively (but optionally), the white matter and CSF time-courses can be extracted from the data and stored as .mat files by putting lines 25 and 26 of SubjectLevelAnalysis.m in comment, than de-commenting lines 27 and 28. Next, enter the following lines of code in the command window:

reslice_mask = 1;
white_matter_parameters = ObtainWhiteMatterPredictor(fileInfo, reslice_mask);
save('white_matter_parameters.mat', 'white_matter_parameters');
CSF_parameters = ObtainCSFPredictor(fileInfo, reslice_mask);
save('CSF_parameters.mat', 'CSF_parameters');

By using this code, the white matter and CSF nuisance parameters will be stored in the white_matter_parameters.mat and CSF_parameters.mat Matlab files, respectively.

Note: the mean white matter and CSF signals are computed from the corresponding segmented anatomical data files, eroded by one voxel in each axis via the ErodeNuisanceVariableMask.m script.

4.2 Subject-level statistical analyses

The subject-level statistical analysis is conducted in a multi-step procedure:
1. The SubjectLevelAnalysis.m script loops over all subjects and stores the GLM statistical beta maps in .nii files in the functional (subjects) data directories.
2. The ROIWiseGLM.m script conducts the actual voxel-wise GLM analysis for the current subject.
3. The ObtainROIPredictor.m script extracts the time-course of the seed- and nuisance ROIs.

4.2.1 Subject-level analyses

Similar to the confirmatory ROI analyses, the subject-level exploratory analyses are governed by the SubjectLevelAnalysis.m script. This function batches the voxel-wise GLM analyses performed by the VoxelWiseBetaMap.m script. It then stores the resulting beta-weight maps of the predictor of interest (i.e., the seed-region) as a .nii file in the subject functional scan directories. 

Optional: It is possible to enter a second brain region as nuisance variable in the voxel-ROI-wise GLM analyses by entering an optional additional input argument (i.e., via the definition statement varargin in Matlab) for the SubjectLevelAnalysis.m and ROIWiseGLM.m functions. This approach is one of two methods to control for partial volume sharing with a region-of-no-interest provided by the analysis pipeline (see section 3.2.3 for the second method; see also Varkevisser, Gladwin, Heesink, van Honk, & Geuze, 2017). 

Optional: As already alluded to in section 3.1, the subject-level GLM analyses are conducted with the (filtered) white matter and CSF parameters as nuisance variables. The SubjectLevelAnalysis.m script can either extract these nuisance variables automatically (if lines 26 and 27 are in comment, but lines 24 and 25 are de-commented) or load these parameters from pre-computed .mat files (if lines 26 and 27 are de-commented, but lines 24 and 25 are in comment). If the former of these options is selected, then the ObtainWhiteMatterPredictor.m and ObtainCSFPredictor.m scripts are executed automatically by the SubjectLevelAnalysis.m function (see also section 3.1).

4.2.1.1 Basolateral amygdala

To execute the subject-level exploratory voxel-wise GLM procedure for the basolateral amygdala as seed-region and the centromedial amygdala as nuisance variable, enter the following lines of code in the command window:

seed_region = 'rAmygdala_LB.nii';
nuisance_region = 'rAmygdala_CM.nii';
non_seed_region = 'Not used';
analysis_type = 1;
SubjectLevelAnalysis(seed_region, non_seed_region, analysis_type, fileInfo, nuisance_region);

4.2.1.2 Centromedial amygdala

To execute the subject-level exploratory voxel-wise GLM procedure for the centromedial amygdala as seed-region and the basolateral amygdala as nuisance variable, enter the following lines of code in the command window:

seed_region = 'rAmygdala_CM.nii';
nuisance_region = 'rAmygdala_LB.nii';
non_seed_region = 'Not used';
analysis_type = 1;
SubjectLevelAnalysis(seed_region, non_seed_region, analysis_type, fileInfo, nuisance_region);

4.2.1.3 Anterior cingulate cortex

To execute the subject-level exploratory voxel-wise GLM procedure for the anterior cingulate cortex (ACC) as seed-region, enter the following lines of code in the command window:

seed_region = 'rCingulum_33.nii';
non_seed_region = 'Not used';
analysis_type = 1;
SubjectLevelAnalysis(seed_region, non_seed_region, analysis_type, fileInfo);

4.2.1.4 Anterior insular cortex

To execute the subject-level exploratory voxel-wise GLM procedure for the anterior insular cortex (AIC) as seed-region, enter the following lines of code in the command window:

seed_region = 'rAIC.nii';
non_seed_region = 'Not used';
analysis_type = 1;
SubjectLevelAnalysis(seed_region, non_seed_region, analysis_type, fileInfo);

4.2.2 Voxel-wise GLM analysis

The voxel-wise GLM analysis is conducted using the VoxelWiseBetaMap.m script, which requires a subject number (corresponding to the line of the current subject in the subjects text file) as input argument. This function first extracts the raw time-course of the (mean) seed-region signal via the ObtainROIPredictor.m script. It then filters the signal of the seed-region predictor of interest, as well as that of the nuisance variables, by using the ezfilt.m sub-function. The time-course of each voxel is then extracted and again filtered with the ezfilt.m script. The following GLM is subsequently conducted for each voxel separately:

Voxel = Intercept + Seed + Global Mean + Motion Parameters + White Matter + CSF + error

To execute the VoxelWiseBetaMap.m function for a single subject by hand, enter the following code in the command window (iSubject = 32 is provided as an example):

seed_region = 'rAmygdala_LB.nii';
iSubject = 32;
load('CSF_parameters.mat');
CSF_parameters = CSF_parameters{iSubject};
load('white_matter_parameters.mat');
white_matter_parameters = white_matter_parameters{iSubject};
[correlation_map_right_ROI, correlation_map_left_ROI, H_Scan] = VoxelWiseBetaMap(seed_region, iSubject, fileInfo, CSF_parameters, white_matter_parameters);

Note: The ezfilt.m script is a function created by Thomas Gladwin (https://www.tegladwin.com). It high- (output argument =  fast) and low-pass (output argument = slow) filters the variable signal at the specified critical value (cutFreq), using a frequency (Fs) defined as 1/TR (i.e., the repetition time of the scans). For the present purposes, the signal will be band-pass filtered at 0.01-0.08 Hz. An example of how ezfilt.m can be used for a given signal is provided below:

TR = 1.6;
[~, signal] = ezfilt(signal, 1/TR, 0.01);   
[signal, ~] = ezfilt(signal, 1/TR, 0.08);

4.2.3 Extraction of the predictors-of-interest

As mentioned section 3.2.3, the raw time-series of the predictor variables of interest are extracted via the ObtainROIPredictor.m script. This function reads the probability/binary atlas maps of the (seed-) regions-of-interest and extracts the BOLD signal of this ROI for each hemisphere separately. 

Optional: The ObtainROIPredictor.m script provides the second option to control for the effects of partial volume sharing via the definition statement varargin (see also section 3.2.1). More specifically, by entering an optional additional nuisance region as input argument, the function is able to control for partial volume sharing by selecting only those voxels for the seed-region, for which the rule seed_probability_map > nuisance_probability_map applies (see also Varkevisser, Gladwin, Heesink, van Honk, & Geuze, 2017). 

To execute the ObtainROIPredictor.m function for a single subject by hand, enter the following code in the command window (iSubject = 19 is provided as an example):

iSubject = 19;
nScans = length(fileInfo.functional_file_names{iSubject});
for iScan = 1:nScans
    H_Scan = spm_vol(fileInfo.functional_file_names{iSubject}{iScan});
    D_Scan = spm_read_vols(H_Scan);
    all_subject_scans{iScan} = D_Scan;
end	
region_of_interest = 'rAmygdala_LB.nii';
nuisance_region = 'rAmygdala_CM.nii';
[predictor_vector_right_ROI, predictor_vector_left_ROI] = ObtainROIPredictor(fileInfo, region_of_interest, all_subject_scans, nuisance_region);

4.3 Group-level statistical analyses

The group-level statistical analyses are governed by a script called GroupLevelAnalysis.m. In this function, the exploratory group ROI analyses are conducted for a given seed-region via voxel-wise independent samples t-tests performed on the subject-level beta-weight maps (see section 4.2). The statistical output is stored as a set of within and between groups t-maps (.nii files). 

4.3.1 Basolateral amygdala

To execute an exploratory group analysis with the basolateral amygdala as seed-region, enter the following code in the command window:

seed_region = 'rAmygdala_LB.nii';
non_seed_region = 'Not used';
analysis_type = 1;
GroupLevelAnalysis(fileInfo, seed_region, non_seed_region, analysis_type);

4.3.2 Centromedial amygdala

To execute an exploratory group analysis with the basolateral amygdala as seed-region, enter the following code in the command window:

seed_region = 'rAmygdala_CM.nii';
non_seed_region = 'Not used';
analysis_type = 1;
GroupLevelAnalysis(fileInfo, seed_region, non_seed_region, analysis_type);

4.3.3 Anterior cingulate cortex

To execute an exploratory group analysis with the anterior cingulate cortex (ACC) as seed-region, enter the following code in the command window:

seed_region = 'rCingulum_33.nii';
non_seed_region = 'Not used';
analysis_type = 1;
GroupLevelAnalysis(fileInfo, seed_region, non_seed_region, analysis_type);

4.3.4 Anterior insular cortex

To execute an exploratory group analysis with the anterior cingulate cortex (ACC) as seed-region, enter the following code in the command window:

seed_region = 'rAIC.nii';
non_seed_region = 'Not used';
analysis_type = 1;
GroupLevelAnalysis(fileInfo, seed_region, non_seed_region, analysis_type);

4.4 Cluster analysis

The GroupLevelAnalysis.m script employs a clustering algorithm in order to localize clusters of voxel-wise significant t-values (at some specified threshold of t, degrees-of-freedom, and minimum voxel extent). This clustering is performed according to the method of Lieberman & Cunningham (2009), which adopts quite liberal criteria for cluster-wise significance. The clustering of the voxel-wise t-maps is governed by a script called LiebermanCunninghamMethod.m, which either runs automatically, i.e., as part of the GroupLevelAnalysis.m function, or can be called by hand.

To execute a cluster analysis of a given voxel-wise within- or between-subjects t-map by hand, enter the following lines of code in the command window (as per example):

input_filename = 'Between_samples_t_map_left_rAmygdala_LB.nii';
input_cluster_criteria.p_thres = 0.0025;
input_cluster_criteria.df = 56;
input_cluster_criteria.min_vox = 20;
output_filename = LiebermanCunninghamMethod(input_filename, input_cluster_criteria);

5. Software testing

The software pipeline contains a data simulation module that enables you to test and validate the analysis scripts. This data simulation module is governed by the DataSimulation.m script, which can either be run automatically via the FileOrganizer.m script (recommended), or by hand (not recommended). 

To run the data simulation automatically via  FileOrganizer.m, enter the following lines of code in the command window:

simulate_data = 1;
fileInfo = FileOrganizer(simulate_data);

To run the data simulation by hand for one particular subject via the DataSimulation.m script, use the following lines of code as an example:

fileInfo = FileOrganizer(0);

iSubject = 25;
this_subject_group = fileInfo.all_subjects{2}{iSubject};
all_functional_scans = fileInfo.functional_file_names{iSubject};
for iScan = 1:length(all_functional_scans)
    this_scan_name = all_functional_scans{iScan};
    
	  simulated_data_dir = [fileInfo.data_dir '\Test simulation'];
	  this_scan_name = DataSimulation(iScan, this_scan_name, 			this_subject_group, fileInfo.base_dir, simulated_data_dir);
end

To test and/or validate the software pipeline, rerun all above-described analysis steps using the newly generated simulated data file information stored in the fileInfo structural variable.
