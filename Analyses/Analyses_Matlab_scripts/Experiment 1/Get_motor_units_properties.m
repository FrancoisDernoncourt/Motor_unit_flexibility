
clc;
close("all")
clear("all")
set(0,'DefaultFigureWindowStyle','normal')
set(0,'DefaultFigureWindowState','maximized')

corresponding_columns_in_MUs_matched = [4,5,6,7];
reference_condition = 1; % which index in "corresponding_columns_in_MUs_matched" is the reference condition
muscle = "GM";
override_subject_idx = true; % if false, automatically picks the 2nd character of the file as subject idx
if override_subject_idx
    subject_idx = "60";
end

MVC_percent = 20;
fsamp = 2048;

% rule for recruitment threshold = first spike of the first series of 3 spikes happening within
% XXX ms
first_spike_burst_max_time = 0.5; % in seconds

nb_of_conditions = numel(corresponding_columns_in_MUs_matched);

files_each_condition = {};
window_lims_files_each_condition = {};

grids_to_load = 1:4; %1:6; %1:4;
grids_to_discard = []; % [1,2];
muscle_per_grid = {"GM";"GM";"GM";"GM"}; %{"RF";"VM";"VL";"VL";"VL";"VL"};
condition_names = {"plateau";"sin0.25";"sin1";"sin3"};

paths_string_to_load = {};
files_string_to_load = {};

fsamp = 2048;

%% MU PROPERTIES - for each trial, and mean over trials
% Recruitment threshold (in % MVC)
% Derecruitment threshold (in % MVC)
% RT rank (relative to other MUs)
% Mean DR ( 1 / mean(ISI) )
% Max DR
% Status (matched / not matched / intermittent)

%% GET PATHS AND FILENAMES

% signal files
for filei=1:nb_of_conditions
    [files_string_to_load{end+1}, paths_string_to_load{end+1}] = ...
        uigetfile(".mat",strcat("Please load the signal file corresponding to condition #",num2str(filei)), ...
        strcat("Please load the signal file corresponding to condition #",num2str(filei)), ...
        "Multiselect","off");
end

if ~override_subject_idx
    subject_idx = files_string_to_load{1};
    subject_idx = subject_idx(2);
end
savefilename = strcat("S",subject_idx,"_",muscle);

% window limitis files
for filei=1:nb_of_conditions
    [files_string_to_load{end+1}, paths_string_to_load{end+1}] = ...
        uigetfile(".mat",strcat("Please load the file with window edges corresponding to condition #",num2str(filei)), ...
        strcat("Please load the signal file window edges corresponding to condition #",num2str(filei)), ...
        "Multiselect","off");
end

% matched file
[files_string_to_load{end+1}, paths_string_to_load{end+1}] = ...
    uigetfile(".mat",strcat("Please load the matched MUs file"), ...
    strcat("Please load the matched MUs file"), ...
    "Multiselect","off");

% FA output file
[files_string_to_load{end+1}, paths_string_to_load{end+1}] = ...
    uigetfile(".mat",strcat("Please load the FA output file"), ...
    strcat("Please load the FA output file"), ...
    "Multiselect","off");

%% ACTUALLY LOAD FILES
load(strcat(paths_string_to_load{end-1},files_string_to_load{end-1})); % Load MUs matched file
load(strcat(paths_string_to_load{end},files_string_to_load{end})); % Load FA output file

grids_to_load_initial = grids_to_load;
grids_to_load = {};

for filei=1:numel(files_string_to_load)
    if filei <= nb_of_conditions % signal files
        load(strcat(paths_string_to_load{filei},files_string_to_load{filei}));
        % Find empty grids (if any)
        grids_to_load{filei} = grids_to_load_initial;
        empty_grids = cellfun(@isempty, edition.Pulsetrain);
        grids_empty_to_discard = find(empty_grids);
        grids_to_discard = [grids_to_discard, grids_empty_to_discard];
        if ~isempty(grids_to_discard)
            for gridi_discard = 1:numel(grids_to_discard)
                grids_to_load{filei}(find(grids_to_load{filei}==grids_to_discard(gridi_discard))) = [];
            end
        end

        %files_each_condition{filei}.force = signal.path;
        files_each_condition{filei}.force = signal.data(max(grids_to_load_initial)*64+1,:);
        if isfield(edition,"time")
            files_each_condition{filei}.time = edition.time;
        else
            files_each_condition{filei}.time = (1:numel(files_each_condition{filei}.force))./fsamp;
        end
        files_each_condition{filei}.binary_spike_trains = edition.Dischargetimes(grids_to_load{filei},:);
        files_each_condition{filei}.condition = MUs_matched.Properties.VariableNames{corresponding_columns_in_MUs_matched(filei)};
    
    elseif (filei > nb_of_conditions) && (filei <= nb_of_conditions*2) % window lims files
        load(strcat(paths_string_to_load{filei},files_string_to_load{filei}));
        window_lims_files_each_condition{filei - nb_of_conditions}.limofEachWindow = limOfEachWindow;
        nb_of_windows = size(limOfEachWindow,1);
        window_lims_files_each_condition{filei - nb_of_conditions}.nb_of_windows = nb_of_windows;
        % concatenated continuous limOfEachWindow
        limofEachWindow_continuous = [];
        for windowi = 1:nb_of_windows
            if windowi == 1
                limofEachWindow_continuous(windowi,1) = 1;
                limofEachWindow_continuous(windowi,2) = limOfEachWindow(windowi,2)-limOfEachWindow(windowi,1);
            else
                limofEachWindow_continuous(windowi,1) = limofEachWindow_continuous(windowi-1,2)+1;
                limofEachWindow_continuous(windowi,2) = limofEachWindow_continuous(windowi,1) + ...
                    (limOfEachWindow(windowi,2)-limOfEachWindow(windowi,1));
            end
        end
        window_lims_files_each_condition{filei - nb_of_conditions}.limofEachWindow_continuous = limofEachWindow_continuous;

    else % MU match file
        break;
    
    end
end

clearvars signal parameters edition filei winowi savename

%% STORE SIGNAL FOR EACH FILE

% change empty cells to zeros in MUs_matched
MUs_matched = MUs_matched(1:max([MUs_matched{:,"MU_idx"}{:}]),:);
for vari = 1:size(MUs_matched,2)
    MUs_matched{find(cellfun(@isempty, MUs_matched{:,vari})),vari} = {0};
end

for conditioni=1:nb_of_conditions
    nb_of_grids = numel(grids_to_load{conditioni});
    binary_matrix = [];
    mu_total = 0;
    corresponding_MUs_for_condition = [];
    for gridi=1:nb_of_grids
        grid_idx = grids_to_load{conditioni}(gridi);
        nonempty_mu_idx_in_grid = find( ~cellfun(@isempty, files_each_condition{conditioni}.binary_spike_trains(gridi,:) ) );
        for mui=1:numel(nonempty_mu_idx_in_grid)
            mu_idx = nonempty_mu_idx_in_grid(mui);
            mu_total = mu_total + 1;
            binary_matrix(mu_total,:) = zeros(1,length(files_each_condition{conditioni}.time));
            binary_matrix(mu_total, files_each_condition{conditioni}.binary_spike_trains{gridi,mu_idx} ) = 1;
            
            % Corresponding MUs
            temp_idx = grid_idx*100 + mu_idx;
            corresponding_MUs_for_condition(mu_total,1) = temp_idx;
            % 1st column is in file
            % 2nd column is MUs_matched
            % row is MU idx when all MUs are together
            corresponding_MUs_matched_row = ...
                find( [MUs_matched{:,corresponding_columns_in_MUs_matched(conditioni)}{:}] == temp_idx);
            if (~isempty(corresponding_MUs_matched_row)) && ...
                    (numel(corresponding_MUs_matched_row) == 1)
                corresponding_MUs_for_condition(mu_total,2) = corresponding_MUs_matched_row;

            elseif (~isempty(corresponding_MUs_matched_row)) && ...
                (numel(corresponding_MUs_matched_row) > 1) % conflict (MU matching with two MUs)
                temp_nb_of_matches_for_conflict = [];
                for case_conflicti=1:numel(corresponding_MUs_matched_row)
                    temp_nb_of_matches_for_conflict(case_conflicti) = ...
                        MUs_matched.MU_found_in_how_many_files{corresponding_MUs_matched_row(case_conflicti)};
                end
                [~,idx_temp] = max(temp_nb_of_matches_for_conflict);
                corresponding_MUs_matched_row = corresponding_MUs_matched_row(idx_temp(1));
                corresponding_MUs_for_condition(mu_total,2) = corresponding_MUs_matched_row;
            else
                corresponding_MUs_for_condition(mu_total,2) = 0;
            end
            corresponding_MUs_for_condition(mu_total,3) = grid_idx;
        end
    end

    binary_matrix_concatenated = [];
    binary_matrix_each_trial = {};
    limofEachWindow = window_lims_files_each_condition{conditioni}.limofEachWindow;
    limofEachWindow_continuous = window_lims_files_each_condition{conditioni}.limofEachWindow_continuous;
    nb_of_windows = size(limofEachWindow,1);
    for windowi=1:nb_of_windows
        current_length = size(binary_matrix_concatenated,2);
        new_length = length(limofEachWindow(windowi,1):limofEachWindow(windowi,2));
        binary_matrix_concatenated(:, ...
           current_length+1:current_length+new_length) = ...
            binary_matrix(:,limofEachWindow(windowi,1):limofEachWindow(windowi,2));
        binary_matrix_each_trial{windowi} = binary_matrix(:,limofEachWindow(windowi,1):limofEachWindow(windowi,2));
    end

    files_each_condition{conditioni}.window_edges = limofEachWindow;
    files_each_condition{conditioni}.binary_matrix_full_file = binary_matrix;
    files_each_condition{conditioni}.binary_matrix_each_trial =  binary_matrix_each_trial;
    files_each_condition{conditioni}.binary_matrix_concatenated = binary_matrix_concatenated;
    files_each_condition{conditioni}.corresponding_matched_MUs = corresponding_MUs_for_condition;
end



%% DETERMINE RT OF EACH MU FOR EACH FILE

%% % First, get force in % of MVC
moving_average_window = 0.2; % in seconds
convolution_kernel = hann(round(moving_average_window*fsamp));
convolution_kernel = convolution_kernel ./ sum(convolution_kernel); % unit area


force_full = files_each_condition{reference_condition}.force;
force_full = conv(force_full,convolution_kernel,'same');

close(figure(1))
figure(1)
plot(force_full)
title(strcat("Select force corresponding to 0% MVC (condition #",num2str(reference_condition),")"));
[input_x,~] = ginput(2);
force_0 = mean( force_full(round(input_x(1)):round(input_x(2)) ));
close(figure(1))

close(figure(1))
figure(1)
plot(force_full)
title(strcat("Select force corresponding to ",num2str(MVC_percent),"% MVC (condition #",num2str(reference_condition),")"));
[input_x,~] = ginput(2);
force_plateau = mean( force_full(round(input_x(1)):round(input_x(2)) ));
close(figure(1))

for conditioni=1:size(files_each_condition,2)
    force_full = files_each_condition{conditioni}.force;
    force_full = force_full - force_0;
    force_full = (force_full ./ (force_plateau - force_0) ) * MVC_percent;
    force_windowed = [];
    limOfEachWindow = files_each_condition{conditioni}.window_edges;
    nb_windows = size(limOfEachWindow,1);
    for windowi=1:nb_windows
        size_of_window = limOfEachWindow(windowi,2) - limOfEachWindow(windowi,1);
        force_windowed(end+1:end+size_of_window+1) = force_full(limOfEachWindow(windowi,1):limOfEachWindow(windowi,2));
    end
    time=(1:length(force_windowed))./fsamp;

    files_each_condition{conditioni}.force_percentMVC = force_full;
    files_each_condition{conditioni}.force_percentMVC_windowed = force_windowed;
end

%% % Then, get recruitment times (for all trials), and recruitment thresholds (for reference condition)
% Rule for Recruitment Time = first spike of the first series of 3 spikes in a row (3
% spikes in < first_spike_burst_max_time)


for conditioni=1:nb_of_conditions

    files_each_condition{conditioni}.corresponding_matched_MUs_table = table();
    files_each_condition{conditioni}.window_to_consider_spikes = {};
    files_each_condition{conditioni}.recruitment_time_for_each_trial = [];
    files_each_condition{conditioni}.RT_per_trial = [];
    files_each_condition{conditioni}.RT_mean_std_rank = [];
    % re-initialize, so that the code cell can be restarted several times

    nb_windows = size(files_each_condition{conditioni}.window_edges,1);
    temp_RT_for_each_trial = nan(size(files_each_condition{conditioni}.binary_matrix_full_file,1),nb_windows);
    binary_spike_trains_total = files_each_condition{conditioni}.binary_matrix_full_file;

    force_signal_for_RT = files_each_condition{conditioni}.force_percentMVC;
    % Divide into trials
    binary_spike_trains_total_each_trial = {};
    force_signal_for_RT_each_trial = {};

    limOfEachWindow_to_use_for_recruitment = files_each_condition{conditioni}.window_edges;
    if conditioni == reference_condition
        % Add the equivalent of 10s after the end of each window lim to
        % include the ramp down part
        limOfEachWindow_to_use_for_recruitment(:,2) = limOfEachWindow_to_use_for_recruitment(:,2) + fsamp*10; % 10'' window after the end of the plateau (and before the beginning of next window)
        limOfEachWindow_to_use_for_recruitment(limOfEachWindow_to_use_for_recruitment>length(force_signal_for_RT)) = length(force_signal_for_RT);
    end

        for triali=1:nb_windows % trial with rampup + plateau + rampdown
            if conditioni == reference_condition
                if triali==1
                    window_for_trial = [1,limOfEachWindow_to_use_for_recruitment(1,2)];
                else
                    window_for_trial = [limOfEachWindow_to_use_for_recruitment(triali-1,2),limOfEachWindow_to_use_for_recruitment(triali,2)];
                end
            else
                window_for_trial = [limOfEachWindow_to_use_for_recruitment(triali,1),...
                limOfEachWindow_to_use_for_recruitment(triali,2)];
            end
            binary_spike_trains_total_each_trial{triali} = binary_spike_trains_total(:, ...
                window_for_trial(1):window_for_trial(2));
            force_signal_for_RT_each_trial{triali} = force_signal_for_RT(window_for_trial(1):window_for_trial(2));
            nb_MUs_temp = size(binary_spike_trains_total_each_trial{triali},1);
            first_spike_per_mu_for_this_trial = [];
            for mui=1:nb_MUs_temp
                spike_idx_temp = find(binary_spike_trains_total_each_trial{triali}(mui,:));
                % if not enough spike, dismiss (less than 4 spikes)
                if numel(spike_idx_temp) < 3
                    temp_RT_for_each_trial(mui,triali) = nan;
                    first_spike_per_mu_for_this_trial(mui) = nan;
                    continue;
                end
                % take the first spike of the first burst of at least 3
                % spikes happening in < first_spike_burst_max_time
                if isempty(spike_idx_temp)
                    temp_RT_for_each_trial(mui,triali) = nan;
                    first_spike_per_mu_for_this_trial(mui) = nan;
                    continue
                end

                diff_spike_idx_temp = diff(spike_idx_temp);
                first_spike_per_mu_for_this_trial(mui) = 0;
                for spiki=1:numel(diff_spike_idx_temp)-2
                    cumulative_ISI_of_three_spikes = sum([diff_spike_idx_temp(spiki), ...
                        diff_spike_idx_temp(spiki+1), ...
                        diff_spike_idx_temp(spiki+2)]);
                    if cumulative_ISI_of_three_spikes < fsamp*first_spike_burst_max_time % if true, recruitment time has been found
                        first_spike_per_mu_for_this_trial(mui) = spike_idx_temp(spiki);
                        break
                    end
                end

                if first_spike_per_mu_for_this_trial(mui) == 0
                    temp_RT_for_each_trial(mui,triali) = nan;
                else
                    temp_RT_for_each_trial(mui,triali) = force_signal_for_RT_each_trial{triali}(first_spike_per_mu_for_this_trial(mui));
                end
            end

            % Visual check
            close(figure((conditioni*10)+triali))
            figure((conditioni*10)+triali)
            % FORCE SIGNAL
            plot(force_signal_for_RT_each_trial{triali});
            hold on
            for mui=1:nb_MUs_temp
                if mui > length(first_spike_per_mu_for_this_trial) % if the MU doesn't discharge during this trial
                    break;
                end
                color_temp = [.5,.5,.5];
                if first_spike_per_mu_for_this_trial(mui) > 0
                    scatter(first_spike_per_mu_for_this_trial(mui),mui,'o','MarkerEdgeColor','red');
                end
                spike_train_plot_temp = binary_spike_trains_total_each_trial{triali}(mui,:);
                spike_train_plot_temp(spike_train_plot_temp==0) = nan;
                spike_train_plot_temp = spike_train_plot_temp .* mui;
                plot(spike_train_plot_temp,'|','Color',color_temp);
            end
            title(strcat("Sanity check for recruitment thresholds - trial #",num2str(triali), ...
                " (condition #",num2str(conditioni),")"));
            saveas(gcf(),strcat(savefilename," - condition #", num2str(conditioni), ...
                " - Recruitment - trial ",num2str(triali),".png"));
            % end of visual check

            % Get recruitment times (for all conditions & trials)
            recruitment_time_for_each_trial = first_spike_per_mu_for_this_trial' + window_for_trial(1);
            files_each_condition{conditioni}.recruitment_time_for_each_trial(:,triali) = recruitment_time_for_each_trial;
            window_to_consider_spikes = [];
            for mui=1:numel(recruitment_time_for_each_trial)
                window_to_consider_spikes(mui,:) = ...
                    [recruitment_time_for_each_trial(mui), limOfEachWindow_to_use_for_recruitment(triali,2)];
            end
            files_each_condition{conditioni}.window_to_consider_spikes{triali} = window_to_consider_spikes;
        end

        % Get recruitment threshold (for reference condition only)
        if conditioni == reference_condition
            mean_std_RT_each_MU = [];
            for mui=1:size(temp_RT_for_each_trial,1)
                mean_std_RT_each_MU(mui,1) = mean(temp_RT_for_each_trial(mui,:),'omitnan');
                mean_std_RT_each_MU(mui,2) = std(temp_RT_for_each_trial(mui,:),'omitnan');
            end

            files_each_condition{conditioni}.RT_per_trial = temp_RT_for_each_trial; % Assign RT to MUs
            files_each_condition{conditioni}.RT_mean_std_rank = mean_std_RT_each_MU; % Assign RT ranks to MUs

            % RANKING BY RT
            files_each_condition{conditioni}.RT_mean_std_rank(:,3) = NaN;
            [~, ranked_mu_indices] = sort(mean_std_RT_each_MU(:,1), 1, 'ascend');
            for ranki=1:length(ranked_mu_indices)
                files_each_condition{conditioni}.RT_mean_std_rank(ranked_mu_indices(ranki),3) = ranki;
            end
            files_each_condition{conditioni}.corresponding_matched_MUs(:,4) = files_each_condition{conditioni}.RT_mean_std_rank(:,3);
            files_each_condition{conditioni}.corresponding_matched_MUs_table = array2table(...
                files_each_condition{conditioni}.corresponding_matched_MUs,'VariableNames',...
                {'In_file_MU_idx','Matched_MU_idx','Grid','RT_rank'});
        else
            files_each_condition{conditioni}.corresponding_matched_MUs_table = array2table(...
                files_each_condition{conditioni}.corresponding_matched_MUs,'VariableNames',...
                {'In_file_MU_idx','Matched_MU_idx','Grid'});
        end % end of "if reference condition"

        files_each_condition{conditioni}.corresponding_matched_MUs_table = addvars(... % add empty variable for muscle
            files_each_condition{conditioni}.corresponding_matched_MUs_table,...
            strings(size(files_each_condition{conditioni}.corresponding_matched_MUs_table,1),1),...
            'After','Grid','NewVariableNames',{'Muscle'});
        for mui=1:size(files_each_condition{conditioni}.corresponding_matched_MUs_table,1)
            files_each_condition{conditioni}.corresponding_matched_MUs_table{mui,'Muscle'} = muscle_per_grid{...
                files_each_condition{conditioni}.corresponding_matched_MUs_table.Grid(mui)};
        end

end % end of "for each condition"

%% GET MEAN DISCHARGE RATES
for conditioni=1:nb_of_conditions
    files_each_condition{conditioni}.mean_DR_per_trial = [];
    files_each_condition{conditioni}.std_DR_per_trial = [];
    files_each_condition{conditioni}.mean_DR_trials_mean = [];
    files_each_condition{conditioni}.std_DR_trials_mean = [];
    for triali=1:numel(files_each_condition{conditioni}.window_to_consider_spikes)
        for mui=1:size(files_each_condition{conditioni}.window_to_consider_spikes{triali},1)
            % Empty by default
            files_each_condition{conditioni}.mean_DR_per_trial(mui,triali) = nan;
            files_each_condition{conditioni}.std_DR_per_trial(mui,triali) = nan;
            if triali==1
                files_each_condition{conditioni}.mean_DR_trials_mean(mui) = nan;
                files_each_condition{conditioni}.std_DR_trials_mean(mui) = nan;
            end

            % nan if didn't get any first spike for the current trial
            if sum(isnan(...
                    files_each_condition{conditioni}.window_to_consider_spikes{triali}(mui,:)...
                    )) > 0
                continue
            end

            % Get discharge rate
            binary_vector_temp = files_each_condition{conditioni}.binary_matrix_full_file(mui,...
                files_each_condition{conditioni}.window_to_consider_spikes{triali}(mui,1) : ...
                files_each_condition{conditioni}.window_to_consider_spikes{triali}(mui,2));
            cumsum_vector_temp = cumsum(binary_vector_temp);
            [~, idx_tmp ] = unique(cumsum_vector_temp);
            IDR_temp = 1./(diff(idx_tmp)./fsamp);
            files_each_condition{conditioni}.mean_DR_per_trial(mui,triali) = mean(IDR_temp);
            files_each_condition{conditioni}.std_DR_per_trial(mui,triali) = std(IDR_temp);
        end
    end
    files_each_condition{conditioni}.mean_DR_trials_mean = ...
        mean(files_each_condition{conditioni}.mean_DR_per_trial, 2, 'omitnan');
    files_each_condition{conditioni}.std_DR_trials_mean = ...
        mean(files_each_condition{conditioni}.std_DR_per_trial, 2, 'omitnan');
end

%% SAVE .csv WITH DISCHARGE RATE PROPERTIES

discharge_rates_per_condition = table();
variable_names = {'Subject','Condition','MU_idx','Muscle','mean_DR','std_DR'};
empty_row = cell(1,numel(variable_names));
discharge_rates_per_condition = cell2table(empty_row,'VariableNames',variable_names);
row_iter = 0;
for conditioni=1:nb_of_conditions
    for mui=1:size(files_each_condition{conditioni}.mean_DR_trials_mean,1)
        row_iter = row_iter + 1;
        discharge_rates_per_condition{row_iter,:} = empty_row;
        discharge_rates_per_condition{row_iter,'Subject'} = {subject_idx};
        discharge_rates_per_condition{row_iter,'Condition'} = condition_names(conditioni);
        discharge_rates_per_condition{row_iter,'MU_idx'} = {...
            files_each_condition{conditioni}.corresponding_matched_MUs_table{mui,"Matched_MU_idx"}};
        discharge_rates_per_condition{row_iter,'Muscle'} = {...
            files_each_condition{conditioni}.corresponding_matched_MUs_table{mui,"Muscle"}};
        discharge_rates_per_condition{row_iter,'mean_DR'} = {...
            files_each_condition{conditioni}.mean_DR_trials_mean(mui)};
        discharge_rates_per_condition{row_iter,'std_DR'} = {...
            files_each_condition{conditioni}.std_DR_trials_mean(mui)};
    end
end

save(strcat("S",subject_idx,"_",muscle,"_table_output_DR.mat"), ...
    'discharge_rates_per_condition');
writetable(discharge_rates_per_condition,...
    strcat("S",subject_idx,"_",muscle,"_table_output_DR.csv"));

%% CONCATENATE BINARY SPIKE TRAINS OF CORRESPONDING MUs
% Only the ones matched in the reference condition

concatenated_file = struct();
concatenated_file.binary_matrix = []; % Will be filled BY RT RANK
concatenated_file.MUs_idx_in_file_to_use = nan(size( ...
        files_each_condition{reference_condition}.corresponding_matched_MUs, 1),...
        nb_of_conditions);
concatenated_file.force = [];
concatenated_file.condition_edges = [];
concatenated_file.trial_edges = [];

for conditioni=1:nb_of_conditions
    % Edges
    concatenated_file.condition_edges(end+1) = length(files_each_condition{conditioni}.force_percentMVC_windowed);
    for triali=1:size(files_each_condition{conditioni}.window_edges,1)
        concatenated_file.trial_edges(end+1) = ...
            files_each_condition{conditioni}.window_edges(triali,2) - files_each_condition{conditioni}.window_edges(triali,1);
    end
    % Force
    concatenated_file.force = [concatenated_file.force, files_each_condition{conditioni}.force_percentMVC_windowed];
    % Discharge matrix
    concatenated_file.binary_matrix( ...
        1:size( files_each_condition{reference_condition}.corresponding_matched_MUs, 1),...
        end+1:end+size(files_each_condition{conditioni}.binary_matrix_concatenated,2)) = nan;
    if conditioni==reference_condition
        concatenated_file.MUs_idx_in_file_to_use(:,conditioni) = 1:size(files_each_condition{reference_condition}.corresponding_matched_MUs,1);
        concatenated_file.MUs_RT_rank = files_each_condition{reference_condition}.corresponding_matched_MUs_table{:,"RT_rank"};
    else
        for mui=1:size(files_each_condition{reference_condition}.corresponding_matched_MUs,1)
            matched_mu_idx_temp = files_each_condition{reference_condition}.corresponding_matched_MUs_table{mui,"Matched_MU_idx"};
            if matched_mu_idx_temp == 0
                continue
            else
                corresponding_mu_idx_in_current_file = find( ...
                    files_each_condition{conditioni}.corresponding_matched_MUs(:,2) ...
                    == matched_mu_idx_temp);
                if isempty(corresponding_mu_idx_in_current_file)
                    continue
                else
                    concatenated_file.MUs_idx_in_file_to_use(mui,conditioni) = corresponding_mu_idx_in_current_file;
                end
            end
        end
    end % end of getting the corresponding MUs
end % end of "for each condition"

corresponding_time_samples = [0,0];
for conditioni=1:nb_of_conditions
    corresponding_time_samples(1) = corresponding_time_samples(2) + 1;
    corresponding_time_samples(2) = corresponding_time_samples(2) + ...
        size(files_each_condition{conditioni}.binary_matrix_concatenated, 2);
    for mui=1:size(concatenated_file.MUs_idx_in_file_to_use,1)
        corresponding_mu_rank = concatenated_file.MUs_RT_rank(mui);
        if isnan(concatenated_file.MUs_idx_in_file_to_use(mui,conditioni))
            continue
        end
        concatenated_file.binary_matrix(corresponding_mu_rank,...
            corresponding_time_samples(1):corresponding_time_samples(2)) = ...
            files_each_condition{conditioni}.binary_matrix_concatenated(...
            concatenated_file.MUs_idx_in_file_to_use(mui,conditioni),:);
    end
end
% Edges
concatenated_file.condition_edges = cumsum(concatenated_file.condition_edges);
concatenated_file.trial_edges = cumsum(concatenated_file.trial_edges);

%% PLOT THE SPIKE TRAINS, FORCE, AND WINDOWS

close(figure(100))
figure(100)

force_subplot_idx = [1,2,3];
spike_trains_subplot_idx = 1:9;
for subplot_idxi = 1:numel(force_subplot_idx)
    spike_trains_subplot_idx(spike_trains_subplot_idx==force_subplot_idx(subplot_idxi)) = [];
end
time = (1:size(concatenated_file.binary_matrix,2))./fsamp;

subplot(3,3,force_subplot_idx)
% Plot force
plot(time, concatenated_file.force,...
    'Color','red','LineWidth',1);
hold on
box off
xlabel("Time (s)")
ylabel("Force (% MVC)")
colorbar("Visible","off")

xlim([0,time(end)])
ylim([0,30])
for triali=1:numel(concatenated_file.trial_edges)-1
    line([concatenated_file.trial_edges(triali),concatenated_file.trial_edges(triali)]./fsamp, ylim(), ...
        'color','black','Linestyle',':','linew',2);
end
for conditioni=1:numel(concatenated_file.condition_edges)-1
    line([concatenated_file.condition_edges(conditioni),concatenated_file.condition_edges(conditioni)]./fsamp, ylim(), ...
        'color','black','Linestyle','-','linew',2.5);
end

subplot(3,3,spike_trains_subplot_idx)
% Plot spike trains
color_scheme = turbo(100);
max_RT = max(files_each_condition{reference_condition}.RT_mean_std_rank(:,1));

background_match_image = ~isnan(concatenated_file.binary_matrix);
background_match_image = imresize(background_match_image,...
    [size(background_match_image, 1), size(background_match_image, 2) * (1/fsamp)]);
background_match_image = background_match_image .* max_RT;
background_image_handle = imagesc(background_match_image,'AlphaData',0.25);
set(gca,'YDir','normal')
colormap("gray")
hold on

for mui=1:size(concatenated_file.binary_matrix, 1)
    spike_train_plot_temp = concatenated_file.binary_matrix(mui,:);

    spike_train_plot_temp(spike_train_plot_temp==0) = nan;
    spike_train_plot_temp = spike_train_plot_temp .* mui;

    corresponding_mu_temp = find(files_each_condition{reference_condition}.RT_mean_std_rank(:,3) == mui);
    temp_mu_RT = files_each_condition{reference_condition}.RT_mean_std_rank(...
        corresponding_mu_temp,...
        1);
    color_idx_temp = max([1, round((temp_mu_RT/max_RT)*100)]); % Can't go below 1
    color_temp = color_scheme(color_idx_temp,:);
    if (files_each_condition{reference_condition}.corresponding_matched_MUs_table{corresponding_mu_temp,'Muscle'} == "VL") ||...
            (files_each_condition{reference_condition}.corresponding_matched_MUs_table{corresponding_mu_temp,'Muscle'} == "GM")
            color_temp = color_scheme(color_idx_temp,:);
    elseif (files_each_condition{reference_condition}.corresponding_matched_MUs_table{corresponding_mu_temp,'Muscle'} == "RF")
        color_temp = [0.8,0,0.8];
    elseif (files_each_condition{reference_condition}.corresponding_matched_MUs_table{corresponding_mu_temp,'Muscle'} == "VM")
        color_temp = [0,0.65,0.65];
    end
    plot(time, spike_train_plot_temp,'|','Color',color_temp);
    hold on
end

xlabel("Time (s)")
ylabel("MU index, organized by RT rank")

for triali=1:numel(concatenated_file.trial_edges)-1
    line([concatenated_file.trial_edges(triali),concatenated_file.trial_edges(triali)]./fsamp, ylim(), ...
        'color','black','Linestyle',':','linew',2);
end
for conditioni=1:numel(concatenated_file.condition_edges)-1
    line([concatenated_file.condition_edges(conditioni),concatenated_file.condition_edges(conditioni)]./fsamp, ylim(), ...
        'color','black','Linestyle','-','linew',2.5);
end
xlim([0,time(end)])
ylim([0.5,size(concatenated_file.binary_matrix, 1)+0.5])
set(background_image_handle, 'XData', ...
    get(background_image_handle, 'XData') - 0.5); % half a pixel shift
box off
cbar = colorbar();
cbar.Limits = [0,max([20,max_RT])];
colormap(cbar,'turbo')
caxis([0,max([20,max_RT])]);
cbar.Label.String = "Recruitment threshold (% MVC)";

sgtitle(strrep(savefilename,"_"," "))
if muscle == "VL"
    subtitle("MUs in magenta are from RF, MUs in dark cyan are from VM")
end

saveas(gcf(),strcat(savefilename,"_MUs_matched_across_conditions_according_to_RT.png"));

%% SAVE .mat OUTPUT OF RECRUITMENT THRESHOLDS AND RANKS
table_output_RT = files_each_condition{reference_condition}.corresponding_matched_MUs_table;
table_output_RT{:,end+1:end+2} = files_each_condition{reference_condition}.RT_mean_std_rank(:,1:2);
table_output_RT.Properties.VariableNames{end-1} = 'RT_mean';
table_output_RT.Properties.VariableNames{end} = 'RT_std';
save(strcat("S",subject_idx,"_",muscle,"_table_output_RT.mat"), ...
    'table_output_RT');
writetable(table_output_RT,...
    strcat("S",subject_idx,"_",muscle,"_table_output_RT.csv"));


%%

% %% GET STATUS (from FA) of MUs
% 
% 
% 
% if reference_condition
%     files_each_condition{1}.RT_mean_std_rank(:,4) = NaN;
%     files_each_condition{1}.corresponding_matched_MUs(:,4) = NaN;
%     %       Necessarily in plateau
%     %       0 if matched with at least one other condition but not included in FA
%     %       1 if matched with at least one other condition but and included in FA
%     %       -1 if unmatched but included in FA
%     %       -2 if unmatched and NOT included in FA
%     for mui=1:size(MUs_clustering_and_correspondance)
%         switch MUs_clustering_and_correspondance.Status{mui}
%             case "matched_discontinuous"
%                 files_each_condition{1}.RT_mean_std_rank(mui,4) = 0;
%                 files_each_condition{1}.corresponding_matched_MUs(mui,4) = 0;
%             case "matched_continuous"
%                 files_each_condition{1}.RT_mean_std_rank(mui,4) = 1;
%                 files_each_condition{1}.corresponding_matched_MUs(mui,4) = 1;
%             case "unmatched_discontinuous"
%                 files_each_condition{1}.RT_mean_std_rank(mui,4) = -2;
%                 files_each_condition{1}.corresponding_matched_MUs(mui,4) = -2;
%             case "unmatched_continuous"
%                 files_each_condition{1}.RT_mean_std_rank(mui,4) = -1;
%                 files_each_condition{1}.corresponding_matched_MUs(mui,4) = -1;
%             otherwise
%                 warning("wrong string somewhere")
%         end
%     end
% else
%     %       -2 if never matched
%     %       -1 if not matched with plateau
%     %       0 if matched with plateau but not included in FA
%     %       1 if matched with plateau and included in FA
%     for mui=1:size(MUs_clustering_and_correspondance)
%         current_matched_MU_idx_temp = MUs_clustering_and_correspondance{mui,2};
%     end
%     error("Not yet implemented for condition other than plateau (reference)");
% end
% 
% % Do a big table / variable called MU properties where there is the RT, the file index, the matched index,
% % the RT rank, the FA status, and all other properties
