function GroupLevelAnalysis(fileInfo, seed_region, non_seed_region, analysis_type)

% Conduct second (group) level analysis: i.e., i.e., resting-state
% connectivity analysis.
%
% Input argument(seed_region): Specify seed-region to base the connectivity
% analysis on, e.g.,
%       seed_region = 'rAmygdala_Total.nii'
% Input argument(non_seed_region): Specify non-seed region-of-interest to
% base the connectivity analysis on, e.g.,
%       non_seed_region  = 'rOFC.nii'
% Input argument(analysis_type): Specify if voxel-wise or ROI-wise analysis
% is to be conducted
%       1. Voxel-wise (exploratory) analysis
%       2. ROI-wise (confirmatory) analysis
% Subfunctions: ExploratoryClusterAnalysis


% ----- Obtain group memberships for all subjects ----- %
group_memberships = cell2mat(fileInfo.all_subjects{2});

% ----- Conduct group analysis ----- %
switch analysis_type
    case 1 % Voxel-wise (exploratory) analysis: compute beta maps
        % ----- Load template functional scan ----- %
        H = spm_vol(fileInfo.functional_file_names{1}{1});
        D = spm_read_vols(H);
        
        % ----- Extract unique group values ----- %
        groupValues = unique(group_memberships);
        nGroups = length(groupValues);
        
        % ----- Loop over hemispheres (within-subjects statistics) ----- %
        for iHemisphere = 1:2
            if iHemisphere == 1
                hemisphere_name = '_left_';
            elseif iHemisphere == 2
                hemisphere_name = '_right_';
            end
            
            % ----- Compute empty statistics matrices ----- %
            n = cell(1, nGroups);
            S = cell(1, nGroups);
            for iGroup = 1:nGroups
                n{iGroup} = length(find(group_memberships == groupValues(iGroup))); % Group sample size(s)
                S{iGroup} = zeros(size(D));
            end
            S2 = S; 
            M = S; V = S; SEM = S; 
            t_within = S;
            
            % ----- Fill S and S2; loop over subjects ----- %
            for iSubject = 1:length(fileInfo.functional_file_names)
                [this_subject_dir, ~, ~] = fileparts(fileInfo.functional_file_names{iSubject}{1});
                iGroup = group_memberships(iSubject) + 1;
                
                % ----- Load beta map for this hemisphere ----- %
                this_beta_map = [this_subject_dir '\Beta_map' hemisphere_name seed_region(1:end-4) '.nii'];
                H_Map = spm_vol(this_beta_map);
                D_Map = spm_read_vols(H_Map);
                
                % ----- Smooth the beta map ----- %
                spm_smooth(D_Map, D_Map, 2); % Consider adjusting the FWHM
                
                % ----- Fill empty sum and  sum of squared beta matrices ----- %
                S{iGroup} = S{iGroup} + D_Map; % Sum of betas map(s)
                S2{iGroup} = S2{iGroup} + (D_Map .^ 2); % Sum of squared betas map(s)
            end
            
            % ----- Fill M, V, SEM, and t_within ----- %
            for iGroup = 1:length(groupValues)
                if iGroup == 1
                    group_name = 'control_group';
                elseif iGroup == 2
                    group_name = 'aggression_group';
                end
                
                % ----- Compute statistical maps ----- %
                M{iGroup} = S{iGroup} ./ n{iGroup}; % Mean of betas map(s)
                V{iGroup} = (S2{iGroup} - n{iGroup} * (M{iGroup} .^ 2)) ./ (n{iGroup} - 1); % Variance of betas map(s)
                SEM{iGroup} = sqrt(V{iGroup}) ./ n{iGroup}; % Standard error of mean betas map(s)
                t_within{iGroup} = M{iGroup} ./ SEM{iGroup}; % One samples t-map(s)
                
                % ----- Save within subjects t-maps to (.nii) file ----- %
                H.fname = ['Within_samples_t_map_' group_name hemisphere_name seed_region(1:end-4) '.nii'];
                H.pinfo(1) = 1;
                H.dt = [16 0];
                spm_write_vol(H, t_within{iGroup});
            end
            
            % ----- Compute between-subjects t-map ----- %
            p1 = ((n{1} - 1) * V{1}) + ((n{2} - 1) * V{2});
            p2 = n{1} + n{2} - 2;
            p3 = (1 / n{1}) + (1 / n{2});
            SE = sqrt((p1 ./ p2) .* p3);
            t_between = (M{2} - M{1}) ./ SE; % independent samples t map
            
            % ----- Write between subjects t-map to (.nii) file ----- %
            H.fname = ['Between_samples_t_map' hemisphere_name seed_region(1:end-4) '.nii'];
            H.pinfo(1) = 1;
            H.dt = [16 0];
            spm_write_vol(H, t_between);
            
            % ----- Lieberman & Cunningham (2009) cluster analysis ----- %
            cluster_criteria.p_thres = 0.0025;
            cluster_criteria.df = p2;
            cluster_criteria.min_vox = 20;
            cluster_filename = LiebermanCunninghamMethod(H.fname, cluster_criteria);
            
            % ----- Extract cluster labels ----- %
            H_Cluster = spm_vol(cluster_filename);
            D_Cluster = spm_read_vols(H_Cluster);
            clusterLabels = unique(D_Cluster);
            clusterLabels = clusterLabels(clusterLabels ~= 0);
            
            % ----- Save cluster statistics to (.csv) file ----- %
            filename = [fileInfo.base_dir '\Results_exploratory_group_analysis_left_' seed_region(1:end-4) '.csv'];
            fileID = fopen(filename, 'w');
            statisticsMatrix{1, 1} = 'Cluster label';
            statisticsMatrix{1, 2} = 'Aggression mean';
            statisticsMatrix{1, 3} = 'Aggression std';
            statisticsMatrix{1, 4} = 'Control mean';
            statisticsMatrix{1, 5} = 'Control std';
            fprintf(fileID, '%s, %s, %s, %s, %s\n', statisticsMatrix{1, 1:end});
                % ----- Loop over clusters ----- %
            for iCluster = 1:length(clusterLabels)
                this_line = iCluster + 1;
                this_cluster_value = clusterLabels(iCluster);
                statisticsMatrix{this_line, 1} = this_cluster_value;
                
                % ----- Obtain voxel indices of this cluster ----- %
                voxel_indices = find(D_Cluster == this_cluster_value);
                
                % ----- Obtain cluster means ----- %
                for iGroup = 1:length(groupValues)
                    this_cluster_mean = nanmean(M{iGroup}(voxel_indices));
                    this_cluster_std = nanstd(M{iGroup}(voxel_indices));
                    
                    % ----- Fill statistics matrix ----- %
                    if groupValues(iGroup) == 1
                        statisticsMatrix{iCluster + 1, 2} = this_cluster_mean;
                        statisticsMatrix{iCluster + 1, 3} = this_cluster_std;
                    elseif groupValues(iGroup) == 0
                        statisticsMatrix{iCluster + 1, 4} = this_cluster_mean;
                        statisticsMatrix{iCluster + 1, 5} = this_cluster_std;
                    end
                end
                
                % ----- Write statistics to file ----- %
                fprintf(fileID, '%f, %f, %f, %f, %f\n', statisticsMatrix{iCluster + 1, 1:end});
                
                % ----- Compute cluster boxplots and save to file ----- %
                this_figure = figure;
                boxplot([M{1}(voxel_indices) M{2}(voxel_indices)], 'Labels', {'Combat control', 'Impulsive aggression'});
                title(['Cluster #' num2str(this_cluster_value)]);
                ylabel('Mean of betas');
                saveas(this_figure, [fileInfo.base_dir '\Boxplot' hemisphere_name seed_region(1:end-4) '_cluster_#' num2str(this_cluster_value)  '.fig']);
                close(this_figure);
            end
            fclose(fileID);
            
        end
        
        
    case 2 % ROI (confirmatory)analysis: compute statistics matrix
        % ----- Obtain all beta values, but seperated by group ----- %
        control_column = 1;        experimental_column = 2;
        right_seed_vs_right_ROI = 1; right_seed_vs_left_ROI = 3;
        left_seed_vs_right_ROI = 2; left_seed_vs_left_ROI = 4;
        all_beta_values{1, control_column} = [];    all_beta_values{1, experimental_column} = [];
        all_beta_values{2, control_column} = [];    all_beta_values{2, experimental_column} = [];
        all_beta_values{3, control_column} = [];    all_beta_values{3, experimental_column} = [];
        all_beta_values{4, control_column} = [];    all_beta_values{4, experimental_column} = [];
        % ----- Loop over all subjects ----- %
        for iSubject = 1:length(fileInfo.functional_file_names)
            [this_subject_dir, ~, ~] = fileparts(fileInfo.functional_file_names{iSubject}{1});
            
            % ----- Load subject seed with non-seed ROI correlations ----- %
            fileContents = load([this_subject_dir '\Output_confirmatory_analysis_' seed_region(1:end-4) '_vs_' non_seed_region(1:end-4) '.mat']);
            fileFields = fields(fileContents);
            betas_seed_vs_ROI = fileContents.(fileFields{1});
            
            % ----- Store all beta values into single cell ----- %
            if group_memberships(iSubject) == 0
                all_beta_values{right_seed_vs_right_ROI, control_column} = [all_beta_values{right_seed_vs_right_ROI, control_column}; betas_seed_vs_ROI{right_seed_vs_right_ROI, 2}];
                all_beta_values{left_seed_vs_right_ROI, control_column} = [all_beta_values{left_seed_vs_right_ROI, control_column}; betas_seed_vs_ROI{left_seed_vs_right_ROI, 2}];
                all_beta_values{right_seed_vs_left_ROI, control_column} = [all_beta_values{right_seed_vs_left_ROI, control_column}; betas_seed_vs_ROI{right_seed_vs_left_ROI, 2}];
                all_beta_values{left_seed_vs_left_ROI, control_column} = [all_beta_values{left_seed_vs_left_ROI, control_column}; betas_seed_vs_ROI{left_seed_vs_left_ROI, 2}];
            elseif group_memberships(iSubject) == 1
                all_beta_values{right_seed_vs_right_ROI, experimental_column} = [all_beta_values{right_seed_vs_right_ROI, experimental_column}; betas_seed_vs_ROI{right_seed_vs_right_ROI, 2}];
                all_beta_values{left_seed_vs_right_ROI, experimental_column} = [all_beta_values{left_seed_vs_right_ROI, experimental_column}; betas_seed_vs_ROI{left_seed_vs_right_ROI, 2}];
                all_beta_values{right_seed_vs_left_ROI, experimental_column} = [all_beta_values{right_seed_vs_left_ROI, experimental_column}; betas_seed_vs_ROI{right_seed_vs_left_ROI, 2}];
                all_beta_values{left_seed_vs_left_ROI, experimental_column} = [all_beta_values{left_seed_vs_left_ROI, experimental_column}; betas_seed_vs_ROI{left_seed_vs_left_ROI, 2}];
            end
        end
        
        % ----- Write statistical output to .csv file ----- %
        filename = [fileInfo.base_dir '\Results_confirmatory_group_analysis_' seed_region(1:end-4) '_vs_' non_seed_region(1:end-4) '.csv'];
        fileID = fopen(filename, 'w');
        statisticsMatrix{1, 1} = 'Contrast';
        statisticsMatrix{1, 2} = 'Hypothesis';
        statisticsMatrix{1, 3} = 'p-value';
        statisticsMatrix{1, 4} = 't-value';
        statisticsMatrix{1, 5} = 'df';
        statisticsMatrix{1, 6} = '95% CI (lower)';
        statisticsMatrix{1, 7} = '95% CI (upper)';
        statisticsMatrix{1, 8} = 'Aggression mean';
        statisticsMatrix{1, 9} = 'Aggression std.';
        statisticsMatrix{1, 10} = 'Control mean';
        statisticsMatrix{1, 11} = 'Control std.';
        fprintf(fileID, '%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n', statisticsMatrix{1, 1:end});
        
        % Left seed vs. left non-seed ROI
        [h, p, ci, stats] = ttest2(all_beta_values{left_seed_vs_left_ROI, control_column}, all_beta_values{left_seed_vs_left_ROI, experimental_column});
        statisticsMatrix{2, 1} = ['Left_' seed_region(1:end-4) '_vs_left_' non_seed_region(1:end-4)];
        statisticsMatrix{2, 2} = h;
        statisticsMatrix{2, 3} = p;
        statisticsMatrix{2, 4} = stats.tstat;
        statisticsMatrix{2, 5} = stats.df;
        statisticsMatrix{2, 6} = ci(1);
        statisticsMatrix{2, 7} = ci(2);
        statisticsMatrix{2, 8} = nanmean(all_beta_values{left_seed_vs_left_ROI, experimental_column});
        statisticsMatrix{2, 9} = nanstd(all_beta_values{left_seed_vs_left_ROI, experimental_column});
        statisticsMatrix{2, 10} = nanmean(all_beta_values{left_seed_vs_left_ROI, control_column});
        statisticsMatrix{2, 11} = nanstd(all_beta_values{left_seed_vs_left_ROI, control_column});
        fprintf(fileID, '%s, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n', statisticsMatrix{2, 1:end});
        
        % Left seed vs. right non-seed ROI
        [h, p, ci, stats] = ttest2(all_beta_values{left_seed_vs_right_ROI, control_column}, all_beta_values{left_seed_vs_right_ROI, experimental_column});
        statisticsMatrix{3, 1} = ['Left_' seed_region(1:end-4) '_vs_right_' non_seed_region(1:end-4)];
        statisticsMatrix{3, 2} = h;
        statisticsMatrix{3, 3} = p;
        statisticsMatrix{3, 4} = stats.tstat;
        statisticsMatrix{3, 5} = stats.df;
        statisticsMatrix{3, 6} = ci(1);
        statisticsMatrix{3, 7} = ci(2);
        statisticsMatrix{3, 8} = nanmean(all_beta_values{left_seed_vs_right_ROI, experimental_column});
        statisticsMatrix{3, 9} = nanstd(all_beta_values{left_seed_vs_right_ROI, experimental_column});
        statisticsMatrix{3, 10} = nanmean(all_beta_values{left_seed_vs_right_ROI, control_column});
        statisticsMatrix{3, 11} = nanstd(all_beta_values{left_seed_vs_right_ROI, control_column});
        fprintf(fileID, '%s, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n', statisticsMatrix{3, 1:end});
        
        % Right seed vs. left non-seed ROI
        [h, p, ci, stats] = ttest2(all_beta_values{right_seed_vs_left_ROI, control_column}, all_beta_values{right_seed_vs_left_ROI, experimental_column});
        statisticsMatrix{4, 1} = ['Right_' seed_region(1:end-4) '_vs_left_' non_seed_region(1:end-4)];
        statisticsMatrix{4, 2} = h;
        statisticsMatrix{4, 3} = p;
        statisticsMatrix{4, 4} = stats.tstat;
        statisticsMatrix{4, 5} = stats.df;
        statisticsMatrix{4, 6} = ci(1);
        statisticsMatrix{4, 7} = ci(2);
        statisticsMatrix{4, 8} = nanmean(all_beta_values{right_seed_vs_left_ROI, experimental_column});
        statisticsMatrix{4, 9} = nanstd(all_beta_values{right_seed_vs_left_ROI, experimental_column});
        statisticsMatrix{4, 10} = nanmean(all_beta_values{right_seed_vs_left_ROI, control_column});
        statisticsMatrix{4, 11} = nanstd(all_beta_values{right_seed_vs_left_ROI, control_column});
        fprintf(fileID, '%s, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n', statisticsMatrix{4, 1:end});
        
        % Right seed vs. right non-seed ROI
        [h, p, ci, stats] = ttest2(all_beta_values{right_seed_vs_right_ROI, control_column}, all_beta_values{right_seed_vs_right_ROI, experimental_column});
        statisticsMatrix{5, 1} = ['Right_' seed_region(1:end-4) '_vs_right_' non_seed_region(1:end-4)];
        statisticsMatrix{5, 2} = h;
        statisticsMatrix{5, 3} = p;
        statisticsMatrix{5, 4} = stats.tstat;
        statisticsMatrix{5, 5} = stats.df;
        statisticsMatrix{5, 6} = ci(1);
        statisticsMatrix{5, 7} = ci(2);
        statisticsMatrix{5, 8} = nanmean(all_beta_values{right_seed_vs_right_ROI, experimental_column});
        statisticsMatrix{5, 9} = nanstd(all_beta_values{right_seed_vs_right_ROI, experimental_column});
        statisticsMatrix{5, 10} = nanmean(all_beta_values{right_seed_vs_right_ROI, control_column});
        statisticsMatrix{5, 11} = nanstd(all_beta_values{right_seed_vs_right_ROI, control_column});
        fprintf(fileID, '%s, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n', statisticsMatrix{5, 1:end});
        fclose(fileID);
end

end
