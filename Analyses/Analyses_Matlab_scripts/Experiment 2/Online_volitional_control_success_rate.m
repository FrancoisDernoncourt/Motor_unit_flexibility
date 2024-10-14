close all
clear all
clc
set(0,'DefaultFigureWindowStyle','normal')
set(0,'DefaultFigureWindowState','maximized')

%% Parameters
subject = "S6";
muscle = "GM";
grid_position = {'Proximal','Distal'};

force_range = [16,24]; % force range to consider
minimum_window_duration = 2; % in s

% Filter discharge rate
% From DelVecchio code neural modules
fsamp = 2048; %set your fsamp;
Wind_s = 0.4;  % hanning window duration % 0.4 for 2.5hz low-pass; 0.2 for 5hz low-pass
HanningW = 2/round(fsamp*Wind_s)*hann(round(fsamp*Wind_s)); %unitary area

% Time spent in quadrant violating size principle
min_meaningful_RT_diff = 0; % in % MVC

plot_signal_within_force_range = false;

%% Get files
[filename_baseline,path_baseline] = ...
    uigetfile(".mat",strcat("Select the baseline file"), ...
    strcat("Select the baseline file"), ...
    "Multiselect","off");

[filename_force,path_force] = ...
    uigetfile(".mat",strcat("Select file with force values"), ...
    strcat("Select file with force values"), ...
    "Multiselect","off");

[filename_RT,path_RT] = ...
    uigetfile(".mat",strcat("Select file with recruitment thresholds"), ...
    strcat("Select file with recruitment thresholds"), ...
    "Multiselect","off");

[filename_tabletrials,path_tabletrials] = ...
    uigetfile(".mat",strcat("Select file with table of trials"), ...
    strcat("Select file with table of trials"), ...
    "Multiselect","off");

[filename_duplicates,path_duplicates] = ...
    uigetfile(".mat",strcat("Select file with duplicate list"), ...
    strcat("Select file with duplicate list"), ...
    "Multiselect","off");

[filenames,path] = ...
    uigetfile(".mat",strcat("Select files with all trials"), ...
    strcat("Select files with all trials"), ...
    "Multiselect","on");

%% LOAD FILES
disp("Loading files");
files = {};
disp(strcat("  Loading smaller files"));

loaded_file = strcat(path_force,filename_force);
load(loaded_file);

loaded_file = strcat(path_RT,filename_RT);
load(loaded_file);
loaded_file = strcat(path_tabletrials,filename_tabletrials);
load(loaded_file);
loaded_file = strcat(path_duplicates,filename_duplicates);
load(loaded_file);

for filei=1:numel(filenames)+1
    disp(strcat("  Loading file ",num2str(filei)," / ", num2str(numel(filenames)+1)));
    if filei==1
        loaded_file = strcat(path_baseline,filename_baseline);
    else
        loaded_file = strcat(path,filenames{filei-1});
    end
    files{filei} = load(loaded_file);
    % Convert force signal into MVC %
    force_percent = files{filei}.signal.path;
    % reverse force if force is negative
    if abs(min(force_percent(fsamp:end-fsamp))) > abs(max(force_percent(fsamp:end-fsamp)))
        force_percent = force_percent .* -1;
    end
    if max(force_percent) < 1.1 % Some files have their "path" files constrained between 0 and 1 for some reason
        % for these files, ask the user to select a window corresponding to
        % 20% of MVC
        figure()
        plot(force_percent)
        mvc_value_temp = mean(force_range);
        title(strcat("File #",num2str(filei), ...
            " has force encoded between 0 and 1. Please select the boundaries for a ", ...
            num2str(mvc_value_temp),"% MVC window"));
        [x_input_pos,~] = ginput(2);
        close(gcf())
        x_input_pos = sort(x_input_pos,'ascend');
        mvc_value_temp = mean(force_percent(x_input_pos(1):x_input_pos(2)))...
            * (1/(mvc_value_temp/100));
        force_percent = force_percent ./ mvc_value_temp;
        force_percent = force_percent .* 100;
        files{filei}.force_percent = force_percent;
    else
        %force_percent = force_percent - force_baseline;
        force_percent = force_percent ./ force_MVC;
        force_percent = force_percent .* 100;
        files{filei}.force_percent = force_percent;
    end

    %% SMOOTHING DRs - directly here to that we can remove the "edition" struct which takes up a lot of RAM
    MUDRs_mat = zeros(1,length(files{filei}.edition.time));
    MUbinary_mat = zeros(1,length(files{filei}.edition.time));
    mu_iter = 0;
    disp(strcat("  Filtering spike trains of file ", num2str(filei)," / ", num2str(numel(filenames)+1)));
    for gridi = 1:size(files{filei}.edition.Pulsetrain,2)
        for mui = 1:size(files{filei}.edition.Pulsetrain{gridi},1)
            mu_iter = mu_iter + 1;
            spike_times_temp = files{filei}.edition.Dischargetimes{gridi,mui};
            MUbinary_mat(mu_iter,:) = zeros(1,length(files{filei}.edition.time));
            MUbinary_mat(mu_iter,spike_times_temp) = 1;
            MUDRs_mat(mu_iter,:) = filtfilt(HanningW,1,MUbinary_mat(mu_iter,:)*fsamp);
            disp(strcat("      Smoothing MU#",num2str(mu_iter)," / ", ...
                num2str(size(recruitment_thresholds,1))));
            if filei==1
                corresponding_MUs(mu_iter,1) = mu_iter;
                corresponding_MUs(mu_iter,2) = gridi*100+mui;
            end
        end
    end
    files{filei}.binary_mat = MUbinary_mat;
    files{filei}.DR_mat = MUDRs_mat;
    files{filei} = rmfield(files{filei},'edition'); % Free up RAM
    files{filei} = rmfield(files{filei},'signal'); % Free up RAM
end

clear("edition","signal","force_percent","spike_times_temp","MUbinary_mat","MUDRs_mat","mu_iter");

%% REMOVE DISCHARGE RATES OF DUPLICATES & OF UNRELIABLE DECOMP
% DUPLICATES

if iscell(duplicates_table{1,1})
    duplicates_table_original = duplicates_table;
    duplicates_table = {};
    for cellx = 1:size(duplicates_table_original,1)
        for celly = 1:size(duplicates_table_original,2)
            duplicates_table{cellx,celly} = duplicates_table_original{cellx,celly}{:};
        end
    end
    duplicates_table = cell2table(duplicates_table, ...
        'VariableNames',duplicates_table_original.Properties.VariableNames);
    clear duplicates_table_original
end

for mui=1:size(duplicates_table,1)
    if ~isnan(duplicates_table{mui,2})
        duplicates_temp = [duplicates_table{mui,1},duplicates_table{mui,2}];
        mu_duplicate_to_remove = duplicates_temp(duplicates_temp~=duplicates_table{mui,3});
        mu_duplicate_to_remove = corresponding_MUs(...
            corresponding_MUs(:,2)==mu_duplicate_to_remove,...
            1);
        mu_duplicate_to_keep = duplicates_temp(duplicates_temp==duplicates_table{mui,3});
        mu_duplicate_to_keep = corresponding_MUs(...
            corresponding_MUs(:,2)==mu_duplicate_to_keep,...
            1);
        for filei=1:numel(files)
            % REPLACE DUPLICATES WITH NANs
            files{filei}.DR_mat(mu_duplicate_to_remove,:) = nan;
            files{filei}.binary_mat(mu_duplicate_to_remove,:) = nan;
           % files{filei}.DR_mat_within_force_range(mu_duplicate_to_remove,:) = nan;
            % REPLACE DUPLICATE TO REMOVE WITH DUPLICATE TO KEEP
%             files{filei}.DR_mat(mu_duplicate_to_remove,:) = files{filei}.DR_mat(mu_duplicate_to_keep,:);
%             files{filei}.binary_mat(mu_duplicate_to_remove,:) = files{filei}.binary_mat(mu_duplicate_to_keep,:);
%             % files{filei}.DR_mat_within_force_range(mu_duplicate_to_remove,:) = files{filei}.DR_mat_within_force_range(mu_duplicate_to_keep,:);
        end
    end
end

% UNRELIABLE DECOMP
for filei=1:numel(files)-1
    % REPLACE UNRELIABLE DECOMP WITH NAN
    mus_unreliable_temp = List_of_MUs_not_decomposed{filei};
    if ~isempty(mus_unreliable_temp)
        for mui=1:numel(mus_unreliable_temp)
            mus_unreliable_temp(mui) = corresponding_MUs(corresponding_MUs(:,2)==mus_unreliable_temp(mui));
        end
        files{filei}.DR_mat(mus_unreliable_temp,:) = nan;
        files{filei}.binary_mat(mus_unreliable_temp,:) = nan;
    end
end

%% Get windows of interest for each file

for filei=1:numel(files)
    files{filei}.chunks_of_continuous_samples = {};

    samples_above_force_min = files{filei}.force_percent>force_range(1);
    samples_above_force_max = files{filei}.force_percent>force_range(2);
    samples_in_force_range = find(samples_above_force_min - samples_above_force_max);
    files{filei}.samples_in_force_range = samples_in_force_range;

    chunks_of_continuous_samples = {nan};
    for sampli=2:numel(samples_in_force_range)
        if (samples_in_force_range(sampli) - samples_in_force_range(sampli-1)) > 1
            chunks_of_continuous_samples{end+1} = [];
        else
            chunks_of_continuous_samples{end}(end+1) = samples_in_force_range(sampli);
        end
    end
    chunks_of_continuous_samples{1}(1) = []; % remove the initial nan value

    long_chunks_idx = find(cellfun(@numel, chunks_of_continuous_samples)...
        > minimum_window_duration*fsamp); % find windows of more than 2 seconds within force range

    nb_windows = numel(long_chunks_idx);
    if (nb_windows < 1) && plot_signal_within_force_range
        figure()
        plot([0]);
        title(strcat("No valid window found in file #",num2str(filei)))
        continue
    end
    individual_windows = nan(nb_windows,2);
    MU_DRs_mat_window_cut = {};
    files{filei}.DR_mat_within_force_range = [];
    for windowi=1:nb_windows
        individual_windows(windowi,:) = [...
            chunks_of_continuous_samples{long_chunks_idx(windowi)}(1),...
            chunks_of_continuous_samples{long_chunks_idx(windowi)}(end)];
        MU_DRs_mat_window_cut{end+1} = files{filei}.DR_mat(:,...
            individual_windows(windowi,1):individual_windows(windowi,2));
        files{filei}.DR_mat_within_force_range = ...
            [files{filei}.DR_mat_within_force_range, MU_DRs_mat_window_cut{end}];
        files{filei}.chunks_of_continuous_samples{windowi} =  [individual_windows(windowi,1):individual_windows(windowi,2)];
    end

    if plot_signal_within_force_range
        figure()
        max_DR = ceil(max(max(files{filei}.DR_mat)));
        for windowi=1:nb_windows
            subplot(2,nb_windows,windowi)
            plot(MU_DRs_mat_window_cut{windowi}')
            xlim([1,size(MU_DRs_mat_window_cut{windowi},2)])
            ylim([0,max_DR])
        end
        subplot(2,nb_windows,nb_windows+1:nb_windows*2)
        plot(files{filei}.DR_mat_within_force_range')
        sgtitle({strcat("Windows of interest for file #",num2str(filei)); ...
            strcat("Conditions = force ranging ", num2str(force_range(1)),"-",num2str(force_range(2)), ...
            "% MVC, during at least ",num2str(minimum_window_duration),"s")});
    end
end

%% Display state-space plot (according to windows of interest) for the different MU pairs and different quadrant targets (or baseline ; no target)
% Also calculate the % time spent in each quadrant
% Take into account MUs not decomposed properly, and give NaN values
% instead (+ say on the graph that not decomposed)

mu_pairs_colors = lines(6);
mu_pairs_colors(3,:) = [];
mu_pairs_colors(1,:) = [];

table_var_names = {...
    'Subject',...
    'Muscle',...
    'MUx',...
    'MUx_position',...
    'MUx_RT', ...
    'MUy',...
    'MUy_position',...
    'MUy_RT', ...
    'Target', ...
    'Time_percent_upper_quadrant', ...
    'Time_percent_lower_quadrant', ...
    'Time_percent_spent_in_target', ...
    'Time_in_force_range'
    };
empty_row = nan(1,numel(table_var_names));
output_cell_array = num2cell(empty_row);

axis_limits = [0,20];
quadrant_angle = cos(pi/3);
nb_of_pairs = size(unique_pairs,1);
targets = {"Baseline (no target)";"Upper quadrant";"Lower quadrant"};
close(figure(999))
figure(999)
row_iter = 0;

for targeti=1:numel(targets)
    for mu_pairi = 1:nb_of_pairs
        pair_to_consider_temp = unique_pairs(mu_pairi,:);
        % Modify unique_pairs in case one MU is a duplicate
        for mui=1:2
             mu_idx_of_pair(mui) = find([duplicates_table{:,1}]==pair_to_consider_temp(mui));
        end
        if sum(isnan([duplicates_table{mu_idx_of_pair,2}]))<2 % one of the MU is a duplicate
            for mui=1:2
                if ~isnan(duplicates_table{mu_idx_of_pair(mui),3})
                    mu_idx_to_keep_temp = duplicates_table{mu_idx_of_pair(mui),3};
                    pair_to_consider_temp(mui) = mu_idx_to_keep_temp;
                end
            end
        end
        unique_pairs(mu_pairi,:) = pair_to_consider_temp;

        idx_of_pair(1) = find(recruitment_thresholds(:,1)==unique_pairs(mu_pairi,1));
        idx_of_pair(2) = find(recruitment_thresholds(:,1)==unique_pairs(mu_pairi,2));
        switch targets{targeti}
            case "Baseline (no target)"
                file_idx_to_consider = 1;
                target_temp = "None";
            case "Upper quadrant"
                upper_quadrants_trials = find(Quadrants_order=="y");
                file_idx_to_consider = upper_quadrants_trials(mu_pairi)+1;
                target_temp = "upper_quadrant";
            case "Lower quadrant"
                lower_quadrants_trials = find(Quadrants_order=="x");
                file_idx_to_consider = lower_quadrants_trials(mu_pairi)+1;
                target_temp = "lower_quadrant";
        end

        %% SET TABLE VALUES

        % Are the MUs from the pair duplicates ?
        % MUx
        mux_RT = (round(recruitment_thresholds(idx_of_pair(1),2),1));
        muy_RT = (round(recruitment_thresholds(idx_of_pair(2),2),1));

        row_iter = row_iter + 1;
        output_cell_array(row_iter,:) = num2cell(empty_row);
        disp(strcat("row_iter = ", num2str(row_iter), ...
            " ; target = ", num2str(targeti), " ; mu_pair = ", num2str(mu_pairi)))

        output_cell_array{row_iter,1} = subject; %     'Subject'
        output_cell_array{row_iter,2} = muscle; %     'Muscle'
        output_cell_array{row_iter,3} = unique_pairs(mu_pairi,1); %     'MUx'
        %     'MUx_position'
        output_cell_array{row_iter,4} = grid_position{floor(unique_pairs(mu_pairi,1)/100)};
        output_cell_array{row_iter,5} = mux_RT; % 'MUx_RT'
        output_cell_array{row_iter,6} = unique_pairs(mu_pairi,2); %     'MUy'
        %     'MUy_position'
        output_cell_array{row_iter,7} = grid_position{floor(unique_pairs(mu_pairi,2)/100)};
        output_cell_array{row_iter,8} = muy_RT; %          'MUy_RT'
        output_cell_array{row_iter,9} = target_temp; %     'Target'

        % Check if both MUs were decomposed properly
        muX_DR_temp = files{file_idx_to_consider}.DR_mat_within_force_range(idx_of_pair(1),:);
        muY_DR_temp = files{file_idx_to_consider}.DR_mat_within_force_range(idx_of_pair(2),:);
        time_in_samples_temp = length(muX_DR_temp);
        decomp_unreliable = [];
        decomp_unreliable = List_of_MUs_not_decomposed{file_idx_to_consider};
        is_unreliable=0;
        if ~isempty(decomp_unreliable)
            for unreliable_mui=1:numel(decomp_unreliable)
                is_unreliable = is_unreliable + ...
                    sum(decomp_unreliable(unreliable_mui) == unique_pairs(mu_pairi,:));
            end
        end
        if is_unreliable
                decomp_reliable_temp = false;
                percent_upper_quadrant_temp = nan;
                percent_lower_quadrant_temp = nan;
                percent_in_target = nan;
        else
            decomp_reliable_temp = true;
            % Calculate % time spent in quadrants
            samples_upper_quadrant_temp = sum(muY_DR_temp > (muX_DR_temp.*2));
            percent_upper_quadrant_temp = (samples_upper_quadrant_temp / length(muX_DR_temp))*100;
            samples_lower_quadrant_temp = sum(muX_DR_temp > (muY_DR_temp.*2));
            percent_lower_quadrant_temp = (samples_lower_quadrant_temp / length(muX_DR_temp))*100;
            if target_temp == "upper_quadrant"
                percent_in_target = percent_upper_quadrant_temp;
            elseif target_temp == "lower_quadrant"
                percent_in_target = percent_lower_quadrant_temp;
            else
                percent_in_target = nan;
            end

        end

       output_cell_array{row_iter,10} = percent_upper_quadrant_temp; %     'Time_percent_upper_quadrant'
       output_cell_array{row_iter,11} = percent_lower_quadrant_temp; %     'Time_percent_lower_quadrant'
       output_cell_array{row_iter,12} = percent_in_target;           %     'Time_percent_in_target'
       output_cell_array{row_iter,13} = time_in_samples_temp/fsamp;  %     ' Time_in_force_range'

        %% PLOT
        subplot(numel(targets), nb_of_pairs, ((targeti-1)*nb_of_pairs)+mu_pairi)

        plot(muX_DR_temp, muY_DR_temp, ...
            'color',[mu_pairs_colors(mu_pairi,:),0.3],'linewidth',2.5);
        hold on
        xlabel(strcat("Grid ",num2str(floor(unique_pairs(mu_pairi,1)/100)),...
            " (",grid_position{floor(unique_pairs(mu_pairi,1)/100)},")", ...
            " MU ", num2str(round(mod(unique_pairs(mu_pairi,1),100))), ...
            " - RT = ", num2str(round(mux_RT,2)), "% MVC"), ...
            "FontSize",9);
       ylabel(strcat("Grid ",num2str(floor(unique_pairs(mu_pairi,2)/100)),...
           " (",grid_position{floor(unique_pairs(mu_pairi,2)/100)},")", ...
            " MU ", num2str(round(mod(unique_pairs(mu_pairi,2),100))), ...
            " - RT = ", num2str(round(muy_RT,2)), "% MVC"), ...
            "FontSize",9);
       xlim(axis_limits);
       ylim(axis_limits);

       % Quadrant patches
       upper_quadrant_patch = patch( ...
            [0,axis_limits(2)*quadrant_angle,0,0],...
            [0,axis_limits(2),axis_limits(2),0],...
            [1,0.75,0]);
        if target_temp == "upper_quadrant"
            upper_quadrant_patch.FaceAlpha = 0.4;
            upper_quadrant_patch.EdgeColor = [1,0.7,0];
            upper_quadrant_patch.EdgeAlpha = 1;
            upper_quadrant_patch.LineWidth = 3;
            upper_quadrant_patch.LineStyle = ":";
        else
            upper_quadrant_patch.FaceAlpha = 0.2;
            upper_quadrant_patch.EdgeAlpha = 0;
        end
        hold on
        lower_quadrant_patch = patch( ...
            [0,axis_limits(2),axis_limits(2),0],...
            [0,0,axis_limits(2)*quadrant_angle,0],...
            [0,0.75,1]);
        if target_temp == "lower_quadrant"
            lower_quadrant_patch.FaceAlpha = 0.4;
            lower_quadrant_patch.EdgeColor = [0,0.8,1];
            lower_quadrant_patch.EdgeAlpha = 1;
            lower_quadrant_patch.LineWidth = 3;
            lower_quadrant_patch.LineStyle = ":";
        else
            lower_quadrant_patch.FaceAlpha = 0.2;
            lower_quadrant_patch.EdgeAlpha = 0;
        end

        if ~decomp_reliable_temp
            unreliable_decomp_patch = patch([axis_limits(1),axis_limits(2),axis_limits(2),axis_limits(1)], ...
                [axis_limits(1),axis_limits(1),axis_limits(2),axis_limits(2)],...
                [.6,.6,.6]);
            unreliable_decomp_patch.FaceAlpha = 0.8;
            unreliable_decomp_patch.EdgeAlpha = 0;
            text(axis_limits(2)/2,axis_limits(2)/2, ...
                strcat("Decomposition is unreliable for MU(s) ",...
                num2str(decomp_unreliable)),... % targeti-2 for appropriate offset
                "HorizontalAlignment","center","VerticalAlignment","middle",...
                "FontSize",11,"FontWeight","bold","Color",[.2,.2,.2]);
        end

       %axis square
       box off
       title({...
           strcat(targets{targeti});...
           strcat("Total time in force range = ", ...
           num2str(round(time_in_samples_temp/fsamp,1)),"s"); ...
           strcat("upper quadrant = ", num2str(round(percent_upper_quadrant_temp)), ...
           "% ; lower quadrant = ", num2str(round(percent_lower_quadrant_temp)), ...
           "% ; target quadrant = ", num2str(round(percent_in_target)),"%")...
           }, ...
           "FontSize",7.5);
    end
end
sgtitle(strcat(subject," - ",muscle," - Conditions = force ranging ", ...
    num2str(force_range(1)),"-",num2str(force_range(2)), ...
    "% MVC, during at least ",num2str(minimum_window_duration),"s"));
saveas(gcf(),strcat("Time_spent_in_target_ONLY_MU_PAIRs_WITH_FEEDBACK_force_range",num2str(force_range(1)), ...
    "-",num2str(force_range(2)),".png"));

table_output = cell2table(output_cell_array,...
     "VariableNames",table_var_names);
save(strcat("Time_spent_in_target_ONLY_MU_PAIRs_WITH_FEEDBACK_force_range",num2str(force_range(1)), ...
    "-",num2str(force_range(2)),".mat"),"table_output");
writetable(table_output,...
    strcat("Time_spent_in_target_ONLY_MU_PAIRs_WITH_FEEDBACK_force_range",num2str(force_range(1)), ...
    "-",num2str(force_range(2)),".csv"));

%% FOR EACH TRIAL, GET SUCCESS RATE, AND % TIME VIOLATING SIZE PRINCIPLE

corresponding_MUs_duplicates_removed = corresponding_MUs;
for mui=1:size(corresponding_MUs_duplicates_removed,1)
    if ~isnan(duplicates_table{mui,2})
        duplicate_temp = duplicates_table{mui,1:2};
        duplicate_idx_to_keep = find(duplicates_temp==duplicates_table{mui,3});
        duplicate_to_remove = duplicate_temp;
        duplicate_to_remove(duplicate_idx_to_keep) = [];
        if corresponding_MUs_duplicates_removed(mui,2)==duplicate_to_remove
            corresponding_MUs_duplicates_removed(mui,1:2) = nan;
        end
    end
end

nb_MUs = size(corresponding_MUs,1);
mu_pair_and_trial_iter = 0;

var_names = {...
    'Subject',...
    'Muscle',...
    'MUx',...
    'MUy',...
    'MU_pair',...
    'MUx_position',...
    'MUy_position',...
    'MUx_RT',...
    'MUy_RT',...
    'Trial',...
    'Target',...
    'Target_MU_to_activate',...
    'Displayed_as_feedback',...
    'Time_percent_success_rate',...
    'Target_violating_size_principle',...
    'Time_percent_violation',...
    'RT_diff_betwwen_MUs'
    };
cell_output = cell(1,numel(var_names));

for triali = 1:numel(files)
    if triali == 1
        target_temp = "baseline";
        pair_displayed_temp = [nan,nan];
    else
        target_temp = Quadrants_order{triali-1};
        pair_displayed_temp = unique_pairs(ceil((triali-1)/2),:);
    end
    disp(strcat("target: ",target_temp,"    ;   pair displayed: ",num2str(pair_displayed_temp)));

    files{triali}.target = target_temp;
    files{triali}.pair_displayed = pair_displayed_temp;

    for muX=1:nb_MUs-1
        for muY = muX+1:nb_MUs
            mu_pair_and_trial_iter = mu_pair_and_trial_iter + 1;
            cell_output(mu_pair_and_trial_iter,:) = cell(1,numel(var_names));
            % Fill-in table
            cell_output{mu_pair_and_trial_iter,1} = subject;
            cell_output{mu_pair_and_trial_iter,2} = muscle;
            cell_output{mu_pair_and_trial_iter,3} = corresponding_MUs_duplicates_removed(muX,1);
            cell_output{mu_pair_and_trial_iter,4} = corresponding_MUs_duplicates_removed(muY,1);
            cell_output{mu_pair_and_trial_iter,5} = corresponding_MUs_duplicates_removed(muX,1)*100+corresponding_MUs_duplicates_removed(muY,1); % will be NaN if duplicate
            if ~isnan(corresponding_MUs_duplicates_removed(muX,1))
                cell_output{mu_pair_and_trial_iter,6} = grid_position{floor(corresponding_MUs_duplicates_removed(muX,2)/100)};
                cell_output{mu_pair_and_trial_iter,8} = recruitment_thresholds(muX,2);
            end
            if ~isnan(corresponding_MUs_duplicates_removed(muY,1))
                cell_output{mu_pair_and_trial_iter,7} = grid_position{floor(corresponding_MUs_duplicates_removed(muY,2)/100)};
                cell_output{mu_pair_and_trial_iter,9} = recruitment_thresholds(muY,2);
            end
            cell_output{mu_pair_and_trial_iter,10} = triali;
            cell_output{mu_pair_and_trial_iter,11} = target_temp;
            % column 12 = 'Target_MU_to_activate'
            cell_output{mu_pair_and_trial_iter,12} = nan;
            corresponding_mu_target_idx = [];
            if sum(isnan(pair_displayed_temp))<2
                grids_of_pair_displayed = grid_position(floor(pair_displayed_temp./100));
                if target_temp == 'y'
                    corresponding_mu_target_idx = find(contains(grids_of_pair_displayed,'Distal'));
                    if ~isempty(corresponding_mu_target_idx)
                        corresponding_mu_target_idx = pair_displayed_temp(corresponding_mu_target_idx);
                        corresponding_mu_target_idx = find(corresponding_MUs_duplicates_removed(:,2)==corresponding_mu_target_idx(1)); %(1) is case of a duplicate, because then 'corresponding_mu_target_idx' can have 2 elements
                        cell_output{mu_pair_and_trial_iter,12} = corresponding_mu_target_idx;
                    end
                elseif target_temp == 'x'
                    corresponding_mu_target_idx = find(contains(grids_of_pair_displayed,'Proximal'));
                    if ~isempty(corresponding_mu_target_idx)
                        corresponding_mu_target_idx = pair_displayed_temp(corresponding_mu_target_idx);
                        corresponding_mu_target_idx = find(corresponding_MUs_duplicates_removed(:,2)==corresponding_mu_target_idx(1)); %(1) is case of a duplicate, because then 'corresponding_mu_target_idx' can have 2 elements
                        cell_output{mu_pair_and_trial_iter,12} = corresponding_mu_target_idx;
                    end
                end
            end
            % Is MU pair displayed as feedback?
            feedback_pair = [nan, nan];
            for mui = 1:2
                temp_idx = find(pair_displayed_temp(mui)==corresponding_MUs_duplicates_removed(:,2));
                if ~isempty(temp_idx)
                    feedback_pair(mui) = temp_idx;
                end
            end
            feedback_pair = unique(feedback_pair); % to prevent issues with duplicates
            if (sum(muX==feedback_pair)+sum(muY==feedback_pair))>=2
                cell_output{mu_pair_and_trial_iter,13} = true;
            else
                cell_output{mu_pair_and_trial_iter,13} = false;
            end
            % Get the firing rates - necessary for the rest of the loop
            total_samples = size(files{triali}.DR_mat_within_force_range,2);
            muX_DR_temp = nan(1,total_samples);
            muY_DR_temp = nan(1,total_samples);
            if ~isnan(corresponding_MUs_duplicates_removed(muX,1))
                muX_DR_temp = files{triali}.DR_mat_within_force_range(muX,:);
            end
            if ~isnan(corresponding_MUs_duplicates_removed(muY,1))
                muY_DR_temp = files{triali}.DR_mat_within_force_range(muY,:);
            end
            % Success rate
            if ~isempty(corresponding_mu_target_idx)
                temp_idx = find(corresponding_mu_target_idx==[...
                    corresponding_MUs_duplicates_removed(muX,1),corresponding_MUs_duplicates_removed(muY,1)]);
                if temp_idx==1 % success comes from muX
                    success_samples = sum(muX_DR_temp > (muY_DR_temp.*2));
                elseif temp_idx==2 % success comes from muY
                    success_samples = sum(muY_DR_temp > (muX_DR_temp.*2));
                else
                    success_samples = nan;
                end
                cell_output{mu_pair_and_trial_iter,14} = (success_samples/total_samples)*100;
            else
                cell_output{mu_pair_and_trial_iter,14} = nan;
            end
            % Which quadrant corresponds to a size principle violation
            mux_RT = cell_output{mu_pair_and_trial_iter,8}; %RT of muX
            muy_RT = cell_output{mu_pair_and_trial_iter,9}; %RT of muY
            violating_quadrant = nan;
            if sum(isnan([mux_RT,muy_RT])) <= 0
                sorted_RTs = sort([mux_RT,muy_RT],'ascend');
                if mux_RT == sorted_RTs(1) % muX is lower threshold
                    violating_quadrant = 'y'; % upper quadrant (distal MU) is the size-principle violating quadrant
                else
                    violating_quadrant = 'x'; % lower quadrant (proximal MU) is the size-principle violating quadrant
                end
            end
            cell_output{mu_pair_and_trial_iter,15} = violating_quadrant;
            % Percent of time spent violating the size principle
            %   % cell_output{mu_pair_and_trial_iter,6} % position (proximal or distal) corresponding to muX
            %   % cell_output{mu_pair_and_trial_iter,7} % position (proximal or distal) corresponding to muY
            if ~isnan(violating_quadrant)
                switch violating_quadrant
                    case 'y'
                        samples_violating_size_principle = sum(muY_DR_temp > (muX_DR_temp.*2));
                    case 'x'
                        samples_violating_size_principle = sum(muX_DR_temp > (muY_DR_temp.*2));
                end
            else
                samples_violating_size_principle = nan;
            end
            cell_output{mu_pair_and_trial_iter,16} = (samples_violating_size_principle/total_samples)*100;
            % Difference of recruitment thresholds
            cell_output{mu_pair_and_trial_iter,17} = abs(mux_RT - muy_RT);

        end % end of "for each mu X"
    end % end of "for each mu Y"
end % end of "for each file"

% convert all empty cells to NaNs
for celli = 1:numel(cell_output)
    if isempty(cell_output{celli})
        cell_output{celli} = nan;
    end
end
table_output = cell2table(cell_output,"VariableNames",var_names);

save(strcat("Table_output_force_range_",num2str(force_range(1)), ...
    "-",num2str(force_range(2)),".mat"),"table_output");
writetable(table_output,...
    strcat("Table_output_force_range_",num2str(force_range(1)), ...
    "-",num2str(force_range(2)),".csv"));
