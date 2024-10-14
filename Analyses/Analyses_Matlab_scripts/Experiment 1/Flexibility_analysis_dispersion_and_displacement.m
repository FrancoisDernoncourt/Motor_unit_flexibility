clc;
close("all")
clear("all")
set(0,'DefaultFigureWindowStyle','docked')

title_prefix = "S7_Quadriceps";
subject_idx = str2double(title_prefix{1}(2));
title_suffix = "";
condition_names = {"plateau";"sin0.25";"sin1";"sin3"}; %{"plateau";"sin0.25";"sin1";"sin3"};
empty_conditions = []; % index of empty condition (if any)
all_conditions_concatenated = true;
corresponding_columns_in_MUs_matched = [4,5,6,7]; %[4,5,6]; %[4,5,6,7]; # If empty_conditions ~= [] = X, then
% the "corresponding columns in MUs matched" of index X will be removed
nb_of_conditions = numel(corresponding_columns_in_MUs_matched);
force_reversed = false;
grids_to_load = 1:6; %1:6; %1:4;
grids_to_discard = []; %[1,2]; %[];
muscles_per_grid = {"RF";"VM";"VL";"VL";"VL";"VL"};

use_only_MUs_included_in_FA = false; % % % % % DO NOT FORGET TO SET
if use_only_MUs_included_in_FA
    title_prefix = strcat(title_prefix,"_only_FA_MUs");
end

compute_flexibility_for_entire_population = false;
compute_flexibility_per_MU_pair = true; % % % % % DO NOT FORGET TO SET % Long to compute
if compute_flexibility_per_MU_pair
    title_prefix = strcat(title_prefix,"_MU_pairs");
else
    title_prefix = strcat(title_prefix,"_MU_pop");
end
perform_displacement_calculation_per_mu_pair = true; % long to compute. Happens only if downsampling=true (otherwise waaaaaay too long)

% Saving options
output_folder = title_prefix;

% Filter discharge rate
% From DelVecchio code neural modules
fsamp = 2048; %set your fsamp;
Wind_s = 0.4;  % hanning window duration % 0.4 for 2.5hz low-pass; 0.2 for 5hz low-pass
HanningW = 2/round(fsamp*Wind_s)*hann(round(fsamp*Wind_s)); %unitary area
detrend_DR_for_each_window = false; % detrend the DR separately for each window
% Basically mostly removes the spike-frequency adaptation

% Flexibility calculation parameters
downsampling = true;
downsampled_by_decimating_or_averaging = false; % true for decimating (keeping every N sample),
% false for averaging (average blocks of N samples)
% (used only if downsampling = true)
downsampling_factor = 100; %100

remove_MUs_with_few_spikes_threshold = 10; % if 0, keeps all MUs
threshold_for_discontinuous_discharge_rate = inf * fsamp; % useful for detrend (if used)
% use "inf" if keeping DR at 0 instead of assigning it "NAN" values

% Display dispersion figures for each MU pair
display_figures_dispersion = false;
save_figures_dispersion = false;

if all_conditions_concatenated
    condition_names{end+1} = "all conditions no plateau";
    condition_names{end+1} = "all conditions";
    added_new_concat_condition = false;
end


%% GET PATHS AND FILENAMES

files_string_to_load = {};
paths_string_to_load = {};

if ~isempty(empty_conditions)
    corresponding_columns_in_MUs_matched(empty_conditions) = [];
end

% signal files
for filei=1:nb_of_conditions
    if sum(filei == empty_conditions) < 1
        [files_string_to_load{end+1}, paths_string_to_load{end+1}] = ...
            uigetfile(".mat",strcat("Please load the signal file corresponding to condition #",num2str(filei)), ...
            strcat("Please load the signal file corresponding to condition #",num2str(filei)), ...
            "Multiselect","off");
    else
        paths_string_to_load{end+1} = [];
        files_string_to_load{end+1} = [];
    end
end

% window limitis files
for filei=1:nb_of_conditions
    if sum(filei == empty_conditions) < 1
        [files_string_to_load{end+1}, paths_string_to_load{end+1}] = ...
            uigetfile(".mat",strcat("Please load the file with window edges corresponding to condition #",num2str(filei)), ...
            strcat("Please load the signal file window edges corresponding to condition #",num2str(filei)), ...
            "Multiselect","off");
    else
        paths_string_to_load{end+1} = [];
        files_string_to_load{end+1} = [];
    end
end

% matched file
[files_string_to_load{end+1}, paths_string_to_load{end+1}] = ...
    uigetfile(".mat",strcat("Please load the matched MUs file"), ...
    strcat("Please load the matched MUs file"), ...
    "Multiselect","off");

% FA output (for clustering correspondance and mu_pair_table)
[files_string_to_load{end+1}, paths_string_to_load{end+1}] = ...
    uigetfile(".mat",strcat("Please load the FA output (for the ref condition)"), ...
    strcat("Please load the FA output (for the ref condition)"), ...
    "Multiselect","off");


%% ACTUALLY LOAD FILES
load(strcat(paths_string_to_load{end-1},files_string_to_load{end-1})); % % Load MUs matched file
load(strcat(paths_string_to_load{end},files_string_to_load{end})); % % Load FA output

% Matched MUs %%%%%%%%%%%%%%%%
% If loaded MUs_matched has more entries/columns than there are conditions
% (corresponds to purposely leaving out conditions), remove the
% corresponding columns
columns_for_conditions_in_MUsMatched_file = min(corresponding_columns_in_MUs_matched):(size(MUs_matched,2)-1); %-1 because last column is the number of matched occurences for each MU
col_rm_temp = setdiff(columns_for_conditions_in_MUsMatched_file, corresponding_columns_in_MUs_matched);
if ~isempty(col_rm_temp)
    MUs_matched(:,col_rm_temp) = [];
    % Recompute "found in how many files"
    for mui=1:sum( ~cellfun(@isempty, MUs_matched{:,"MU_idx"})) % discounts
        MUs_matched{mui,"MU_found_in_how_many_files"}{:} = ...
            sum ( ~cellfun(@isempty, MUs_matched{mui,corresponding_columns_in_MUs_matched}) );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(grids_to_discard)
    for gridi_discard = 1:numel(grids_to_discard)
        grids_to_load(find(grids_to_load==grids_to_discard(gridi_discard))) = [];
    end
end

grids_to_load_for_file = {};
for filei=1:numel(files_string_to_load)
    if filei <= nb_of_conditions % signal files

        if isempty(paths_string_to_load{filei})
            % discount condtion if it corresponds to a condition for
            % which there is no file or recording
            files_each_condition{filei}.force = [];
            files_each_condition{filei}.time = [];
            files_each_condition{filei}.binary_spike_trains = [];
            files_each_condition{filei}.condition = condition_names{filei};
            continue;
        end

        load(strcat(paths_string_to_load{filei},files_string_to_load{filei}));
        %files_each_condition{filei}.force = signal.path;

        % Remove empty grids from grids to load
        grids_to_load_for_file{filei} = grids_to_load;
        for gridi=1:numel(grids_to_load)
            if isempty(edition.Pulsetrain{grids_to_load(gridi)})
                grids_to_load_for_file{filei}(find(grids_to_load_for_file{filei}==grids_to_load(gridi))) = [];
            end
        end

        files_each_condition{filei}.force = signal.data(signal.ngrid*64+1,:);
        if force_reversed
            files_each_condition{filei}.force = files_each_condition{filei}.force .* (-1);
        end
        files_each_condition{filei}.time = (1:length(files_each_condition{filei}.force))./fsamp;
        files_each_condition{filei}.binary_spike_trains = edition.Dischargetimes(grids_to_load_for_file{filei},:);
        files_each_condition{filei}.condition = MUs_matched.Properties.VariableNames{corresponding_columns_in_MUs_matched(filei)};

    elseif (filei > nb_of_conditions) && (filei <= nb_of_conditions*2) % window lims files

        if isempty(paths_string_to_load{filei})
            % discount condtion if it corresponds to a condition for
            % which there is no file or recording
            window_lims_files_each_condition{filei - nb_of_conditions}.limofEachWindow = [];
            window_lims_files_each_condition{filei - nb_of_conditions}.nb_of_windows = [];
            window_lims_files_each_condition{filei - nb_of_conditions}.limofEachWindow_continuous = [];
            continue;
        end

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

%% NEED TO UPDATE EMPTY CONDITIONS IF USING ONLY MUS INCLUDED IN FA

if use_only_MUs_included_in_FA
    % get list of MUs to remove
    only_continuous_FA_MUs = MUs_clustering_and_correspondance;
    idx_to_remove_temp = find([only_continuous_FA_MUs.Status{:}] ~= 'matched_continuous');
    only_continuous_FA_MUs(idx_to_remove_temp,:) = [];
    MUs_to_keep_because_present_in_FA = {};

    % find which files will be empty when removing the MUs not used in FA
    for filei=1:numel(corresponding_columns_in_MUs_matched)
        col_temp = corresponding_columns_in_MUs_matched(filei);
        MUs_matched_from_file_temp = []; % these correspond to the MATCHED INDEXES and not the FILE INDEXES
        for mui=1:size(MUs_matched,1)
            if ~isempty(MUs_matched{mui,col_temp}{:})
                MUs_matched_from_file_temp(end+1) = mui;
            end
        end
        idx_matched_to_remove_temp = [];
        for mui=1:numel(MUs_matched_from_file_temp)
            if sum(MUs_matched_from_file_temp(mui)==only_continuous_FA_MUs.MU_idx_in_matched_MUs)<1
                idx_matched_to_remove_temp(end+1) = mui;
            end
        end
        MUs_matched_from_file_temp(idx_matched_to_remove_temp) = [];

        MUs_to_keep_because_present_in_FA{filei} = MUs_matched_from_file_temp;
        % these correspond to the MATCHED INDEXES and not the FILE INDEXES

        if isempty(MUs_matched_from_file_temp) % if no MU remaining, add condition as being empty
            if (isempty(empty_conditions)) || (sum(empty_conditions==filei)==0)
                empty_conditions(end+1) = filei;
            end
        end
    end
end

%% UPDATE MUs MATCHED IF ADDING A CONDITION WITH ALL CONDITIONS CONCATENATED

if all_conditions_concatenated
    new_col_idx = [max(corresponding_columns_in_MUs_matched)+1 , max(corresponding_columns_in_MUs_matched)+2];
    % 2 new cols ; one for all conditions concatenated ; one for all
    % conditions concatenated without plateau
    MUs_matched(:,end+1) = MUs_matched(:,end); % duplicate last column
    MUs_matched(:,end+1) = MUs_matched(:,end); % duplicate last column... twice
    MUs_matched.Properties.VariableNames(end-2) = {condition_names{end-1}{:}};
    MUs_matched.Properties.VariableNames(end-1) = {condition_names{end}{:}};
    MUs_idx_found_in_all_files = [];

    for rowi=1:size(MUs_matched,1)
        if (MUs_matched{rowi,end-1}{:} == (nb_of_conditions - numel(empty_conditions)) ) &&...
                (isempty(intersect(MUs_matched{rowi,"Grid"}{:},grids_to_discard)))
            % this MU is represented in all conditions AND it doesn't belong to
            % a discarded grid
            MUs_matched{rowi,end-1}{:} = [];
            MUs_matched{rowi,end-2}{:} = [];
            for conditioni=1:numel(corresponding_columns_in_MUs_matched)
                MU_to_add_temp = MUs_matched{rowi,corresponding_columns_in_MUs_matched(conditioni)}{:};
                if (sum(conditioni == empty_conditions) >= 1) || isempty(MU_to_add_temp)
                    % Case of an unrepresented condition
                    continue;
                end
                MUs_matched{rowi,end-1}{:}(conditioni) = MUs_matched{rowi,corresponding_columns_in_MUs_matched(conditioni)}{:};
                MUs_matched{rowi,end-2}{:}(conditioni) = MUs_matched{rowi,corresponding_columns_in_MUs_matched(conditioni)}{:};
            end
            MUs_idx_found_in_all_files(end+1) = rowi;
        else
            MUs_matched{rowi,end-1}{:} = [];
            MUs_matched{rowi,end-2}{:} = [];
        end
    end
end

%% CREATE BINARY MATRIX AND SMOOTHED DRs FOR EACH FILE

% change empty cells to zeros in MUs_matched
MUs_matched = MUs_matched(1:max([MUs_matched{:,"MU_idx"}{:}]),:);
for vari = 1:size(MUs_matched,2)
    MUs_matched{find(cellfun(@isempty, MUs_matched{:,vari})),vari} = {0};
end

for conditioni=1:nb_of_conditions

    if sum(conditioni == empty_conditions) >= 1
        % Case of an unrepresented condition
        files_each_condition{conditioni}.binary_matrix = [];
        files_each_condition{conditioni}.binary_matrix_concatenated = [];
        files_each_condition{conditioni}.corresponding_matched_MUs = [];
        files_each_condition{conditioni}.smoothed_DR = [];
        files_each_condition{conditioni}.smoothed_DR_discontinuous_nans = [];
        files_each_condition{conditioni}.smoothed_DR_concatenated_with_nans = [];
        files_each_condition{conditioni}.smoothed_DR_concatenated_no_nans = [];
        files_each_condition{conditioni}.smoothed_DR_per_trial_no_nans = [];
        files_each_condition{conditioni}.smoothed_DR_per_trial_with_nans = [];
        files_each_condition{conditioni}.smoothed_DR_per_trial_downsampled = [];
        continue;
    end

    binary_matrix = [];
    smoothed_DR = [];
    smoothed_DR_discontinuous_nans = [];
    mu_total = 0;
    corresponding_MUs_for_condition = [];
    for gridi=1:numel(grids_to_load_for_file{conditioni})
        grid_idx = grids_to_load_for_file{conditioni}(gridi);
        nonempty_mu_idx_in_grid = find( ~cellfun(@isempty, files_each_condition{conditioni}.binary_spike_trains(gridi,:) ) );
        % At least N spikes are necessary (to be chosen in parameters)
        mu_too_few_spikes = [];
        for mui=1:numel(nonempty_mu_idx_in_grid)
            if numel(files_each_condition{conditioni}.binary_spike_trains{gridi,mui}) <= remove_MUs_with_few_spikes_threshold
                mu_too_few_spikes(end+1) = mui;
            end
        end
        nonempty_mu_idx_in_grid(mu_too_few_spikes) = [];
        clearvars mu_too_few_spikes
        for mui=1:numel(nonempty_mu_idx_in_grid)
            mu_idx = nonempty_mu_idx_in_grid(mui);
            mu_total = mu_total + 1;
            binary_matrix(mu_total,:) = zeros(1,length(files_each_condition{conditioni}.time));
            binary_matrix(mu_total, files_each_condition{conditioni}.binary_spike_trains{gridi,mu_idx} ) = 1;
            smoothed_DR(mu_total,:) = filtfilt(HanningW,1,binary_matrix(mu_total,:)*fsamp);
            % Add NANs to discontinuous gaps
            smoothed_DR_discontinuous_nans(mu_total,:) = smoothed_DR(mu_total,:);
            spike_train_idx = find(binary_matrix(mu_total,:));
            spike_train_diff_temp = diff([0, spike_train_idx, size(binary_matrix,2)]);
            for spiki=1:length(spike_train_diff_temp)
                if spike_train_diff_temp(spiki) > threshold_for_discontinuous_discharge_rate
                    if spiki == 1 % if first spike
                        smoothed_DR_discontinuous_nans(mu_total,1:spike_train_idx(spiki)) = nan;
                    elseif spiki == length(spike_train_diff_temp) % if last spike
                        smoothed_DR_discontinuous_nans(mu_total,spike_train_idx(spiki-1):end) = nan;
                    else
                        smoothed_DR_discontinuous_nans(mu_total,spike_train_idx(spiki-1):spike_train_idx(spiki)) = nan;
                    end
                end
            end
            %

            % Corresponding MUs
            temp_idx = grid_idx*100 + mu_idx;
            corresponding_MUs_for_condition(mu_total,1) = temp_idx;
            % 1st column is in file
            % 2nd column is MUs_matched
            % row is MU idx when all MUs are together
            corresponding_MUs_matched_row = ...
                find( [MUs_matched{:,corresponding_columns_in_MUs_matched(conditioni)}{:}] == temp_idx);
            if ~isempty(corresponding_MUs_matched_row)
                if numel(corresponding_MUs_matched_row) > 1 % In a few files the same MU can be detected twice
                    corresponding_MUs_matched_row = corresponding_MUs_matched_row(1);
                end
                corresponding_MUs_for_condition(mu_total,2) = corresponding_MUs_matched_row;
            else
                corresponding_MUs_for_condition(mu_total,2) = 0;
            end
        end
    end

    binary_matrix_concatenated = [];
    smoothed_DR_concatenated_with_nans = [];
    smoothed_DR_concatenated_no_nans = [];

    limofEachWindow = window_lims_files_each_condition{conditioni}.limofEachWindow;
    limofEachWindow_continuous = window_lims_files_each_condition{conditioni}.limofEachWindow_continuous;
    nb_of_windows = size(limofEachWindow,1);
    for windowi=1:nb_of_windows
        current_length = size(binary_matrix_concatenated,2);
        new_length = length(limofEachWindow(windowi,1):limofEachWindow(windowi,2));
        temp_windowed_binary_mat = binary_matrix(:,limofEachWindow(windowi,1):limofEachWindow(windowi,2));
        binary_matrix_concatenated(:, ...
            current_length+1:current_length+new_length) = ...
            temp_windowed_binary_mat;
        % get inter-spike-interval (ISI) for each window
        % WITH NAN VALUES
        smoothed_DR_current_window = smoothed_DR_discontinuous_nans(:,limofEachWindow(windowi,1):limofEachWindow(windowi,2));
        if detrend_DR_for_each_window
            for mui = 1:size(smoothed_DR_current_window,1)
                smoothed_DR_current_window(mui,:) = ...
                    detrend(smoothed_DR_current_window(mui,:),'omitnan');
            end
        end
        smoothed_DR_concatenated_with_nans(:, ...
            current_length+1:current_length+new_length) = ...
            smoothed_DR_current_window;
        % NO NAN VALUES
        smoothed_DR_current_window = smoothed_DR(:,limofEachWindow(windowi,1):limofEachWindow(windowi,2));
        smoothed_DR_current_window(smoothed_DR_current_window==0) = nan;
        if detrend_DR_for_each_window
            for mui = 1:size(smoothed_DR_current_window,1)
                smoothed_DR_current_window(mui,:) = ...
                    detrend(smoothed_DR_current_window(mui,:),'omitnan');
            end
        end
        smoothed_DR_current_window(isnan(smoothed_DR_current_window)) = 0;
        smoothed_DR_concatenated_no_nans(:, ...
            current_length+1:current_length+new_length) = ...
            smoothed_DR_current_window;
    end
    clearvars temp_max temp_std

    files_each_condition{conditioni}.window_edges = limofEachWindow;
    files_each_condition{conditioni}.window_edges_continuous = limofEachWindow_continuous;
    files_each_condition{conditioni}.binary_matrix = binary_matrix;
    files_each_condition{conditioni}.binary_matrix_concatenated = binary_matrix_concatenated;
    files_each_condition{conditioni}.corresponding_matched_MUs = corresponding_MUs_for_condition;
    % Check the MUs which were used for the FA
    files_each_condition{conditioni}.corresponding_matched_MUs(:,3) = -2;
    %       -2 if never matched
    %       -1 if not matched with plateau
    %       0 if matched with plateau but not included in FA
    %       1 if matched with plateau and included in FA
    for mui = 1:size(corresponding_MUs_for_condition,1)
        matched_idx_temp = corresponding_MUs_for_condition(mui,2);
        if matched_idx_temp > 0
            FA_output_mu_idx_correspondance = find(MUs_clustering_and_correspondance.MU_idx_in_matched_MUs == matched_idx_temp);
            if isempty(FA_output_mu_idx_correspondance)
                % The MU was matched at least once, but not found in the plateau
                % (reference)
                files_each_condition{conditioni}.corresponding_matched_MUs(mui,3) = -1;
            else
                if MUs_clustering_and_correspondance.Status{FA_output_mu_idx_correspondance} == "matched_discontinuous"
                    % The MU was matched with the plateau, but the MU wasn't
                    % used for the FA
                    files_each_condition{conditioni}.corresponding_matched_MUs(mui,3) = 0;
                elseif MUs_clustering_and_correspondance.Status{FA_output_mu_idx_correspondance} == "matched_continuous"
                    % The MU was matched with the plateau, and the MU was used
                    % for the FA
                    files_each_condition{conditioni}.corresponding_matched_MUs(mui,3) = 1;
                end
            end
        end
    end
    %
    files_each_condition{conditioni}.smoothed_DR = smoothed_DR;
    files_each_condition{conditioni}.smoothed_DR_discontinuous_nans = smoothed_DR_discontinuous_nans;
    files_each_condition{conditioni}.smoothed_DR_concatenated_with_nans = smoothed_DR_concatenated_with_nans;
    files_each_condition{conditioni}.smoothed_DR_concatenated_no_nans = smoothed_DR_concatenated_no_nans;

    % Get DR for each trial
    files_each_condition{conditioni}.smoothed_DR_per_trial_no_nans = {};
    files_each_condition{conditioni}.smoothed_DR_per_trial_with_nans = {};
    for triali = 1:size(limofEachWindow_continuous,1)
        files_each_condition{conditioni}.smoothed_DR_per_trial_no_nans{triali} = ...
            smoothed_DR_concatenated_no_nans(:,...
            limofEachWindow_continuous(triali,1):limofEachWindow_continuous(triali,2));
        files_each_condition{conditioni}.smoothed_DR_per_trial_with_nans{triali} = ...
            smoothed_DR_concatenated_with_nans(:,...
            limofEachWindow_continuous(triali,1):limofEachWindow_continuous(triali,2));
    end
end % end of "for each condition"

%% DOWNSAMPLED SIGNAL FOR EACH FILE

if downsampling
    for conditioni=1:nb_of_conditions

        if sum(conditioni == empty_conditions) >= 1
            % Case of an unrepresented condition
            continue;
        end

        files_each_condition{conditioni}.smoothed_DR_per_trial_downsampled = {};
        smoothed_DR_downsampled = [];
        for mui=1:size(files_each_condition{conditioni}.smoothed_DR_concatenated_no_nans,1)
            temp_MU = files_each_condition{conditioni}.smoothed_DR_concatenated_no_nans(mui,:);
            if downsampled_by_decimating_or_averaging
                smoothed_DR_downsampled(mui,:) = downsample(temp_MU, downsampling_factor);
            else
                nb_chunks = floor(length(temp_MU)/downsampling_factor);
                first_sample_chunk = 0;
                last_sample_chunk = 0;
                for sample_chunki=1:(nb_chunks-1)
                    first_sample_chunk = last_sample_chunk + 1;
                    last_sample_chunk = sample_chunki * downsampling_factor;
                    smoothed_DR_downsampled(mui,sample_chunki) = ...
                        mean(temp_MU(first_sample_chunk:last_sample_chunk),'omitnan');
                end
            end
        end
        files_each_condition{conditioni}.smoothed_DR_downsampled = smoothed_DR_downsampled;

        window_lims_for_downsampled_temp = files_each_condition{conditioni}.window_edges_continuous;
        window_lims_for_downsampled_temp(:,1) = ceil(window_lims_for_downsampled_temp(:,1) ./ downsampling_factor) + 1;
        window_lims_for_downsampled_temp(:,2) = floor(window_lims_for_downsampled_temp(:,2) ./ downsampling_factor) - 1;
        window_lims_for_downsampled_temp(end) = size(smoothed_DR_downsampled,2);
        files_each_condition{conditioni}.window_edges_downsampled = window_lims_for_downsampled_temp;
        %"+1 and -1" to be sure that one sample from
        % another trial is not included in the current trial
        for triali=1:size(window_lims_for_downsampled_temp,1)
            files_each_condition{conditioni}.smoothed_DR_per_trial_downsampled{triali} = ...
                smoothed_DR_downsampled(:,...
                (window_lims_for_downsampled_temp(triali,1)):(window_lims_for_downsampled_temp(triali,2)));
        end
    end
end

%% REMOVE MUs NOT USED FOR FA (if necessary, according to selection in "use_only_MUs_included_in_FA" variable)

if use_only_MUs_included_in_FA

    % MUs_to_keep_because_present_in_FA % of note; the indexes are based in
    % the matched MUs index, so they have to be converted into the
    % unmatched indexes

    for conditioni=1:nb_of_conditions
        if sum(empty_conditions==conditioni) > 0
            continue;
        end
        % Remove MUs from each variable
        nb_of_MUs_before_removal = size(files_each_condition{conditioni}.smoothed_DR,1);
        idx_of_no_FA_MUs_for_removal = [];
        MUs_idx_to_keep_current_condition_indexes = ...
            files_each_condition{conditioni}.corresponding_matched_MUs(:,2);
        for mui=nb_of_MUs_before_removal:-1:1
            if sum(MUs_idx_to_keep_current_condition_indexes(mui) == ...
                    MUs_to_keep_because_present_in_FA{conditioni})...
                    < 1
                idx_of_no_FA_MUs_for_removal(end+1) = mui;
            end
        end

        fields_to_check = fieldnames(files_each_condition{conditioni});
        for fieldi = 1:numel(fields_to_check)
            current_field_temp = fields_to_check{fieldi};
            if size(files_each_condition{conditioni}.(current_field_temp),1) == nb_of_MUs_before_removal
                files_each_condition{conditioni}.(current_field_temp)...
                    (idx_of_no_FA_MUs_for_removal,:) = [];
                % something to change there
            elseif contains(current_field_temp,'_per_trial') && iscell(files_each_condition{conditioni}.(current_field_temp))
                % the "for each trial" field
                for triali=1:numel(files_each_condition{conditioni}.(current_field_temp))
                    files_each_condition{conditioni}.(current_field_temp){triali}(idx_of_no_FA_MUs_for_removal,:) = [];
                end
            end
        end % end "for each field"

    end % end "for each condition"

end % end "if use_only_MUs_included_in_FA"

%% CREATE FOLDER
mkdir(output_folder);
for conditioni=1:nb_of_conditions
    mkdir(strcat(output_folder,"/",condition_names{conditioni}));
end
if all_conditions_concatenated
    mkdir(strcat(output_folder,"/",condition_names{end}));
    mkdir(strcat(output_folder,"/",condition_names{end-1}));
end

%% ADD THE CONCATENATED CONDITION, IF NECESSARY
if all_conditions_concatenated && numel(MUs_idx_found_in_all_files)>1

    if ~added_new_concat_condition
        nb_of_conditions = nb_of_conditions+2;
        corresponding_columns_in_MUs_matched(end+1) = max(corresponding_columns_in_MUs_matched)+1;
        corresponding_columns_in_MUs_matched(end+1) = max(corresponding_columns_in_MUs_matched)+1; % do it twice
        added_new_concat_condition = true;
    end
    files_each_condition{nb_of_conditions-1} = struct();
    files_each_condition{nb_of_conditions} = struct();

    % Get list of corresponding MUs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     nb_MUs_matched_in_all_files = sum(...
    %         cellfun(@(x) sum(x)~=0, MUs_matched{:,corresponding_columns_in_MUs_matched(end)})...
    %         );
    nb_MUs_matched_in_all_files = numel(MUs_idx_found_in_all_files); % same as line above, but more simple
    corresponding_MUs_to_fetch = []; % one line per MU, 1st column is corresponding MU in MUs_matched,
    % other columns are corresponding MUs in the different conditions
    corresponding_MUs_to_fetch(:,1) = MUs_idx_found_in_all_files;

    for conditioni=1:(nb_of_conditions-2)

        if sum(conditioni == empty_conditions) >= 1
            % Case of an unrepresented condition
            continue;
        end

        MUs_matched_in_all_conditions_found_in_this_file = [];
        for mui = 1:nb_MUs_matched_in_all_files
            temp_mu_idx = find(...
                files_each_condition{conditioni}.corresponding_matched_MUs(:,2) == MUs_idx_found_in_all_files(mui) );
            if ~isempty(temp_mu_idx)
                MUs_matched_in_all_conditions_found_in_this_file(end+1) = temp_mu_idx;
            else
                MUs_matched_in_all_conditions_found_in_this_file(end+1) = 0;
            end
        end
        corresponding_MUs_to_fetch(:,conditioni+1) = MUs_matched_in_all_conditions_found_in_this_file;
    end

    % Find rows containing at least one zero
    rows_with_zero = any(corresponding_MUs_to_fetch == 0, 2);
    % Remove rows containing zeros
    corresponding_MUs_to_fetch = corresponding_MUs_to_fetch(~rows_with_zero, :);

    colnames = {};
    colnames{1} = 'Matched_MU_idx';
    for coli=2:(nb_of_conditions-1)
        if sum(coli == (empty_conditions+1)) >= 1
            % Case of an unrepresented condition
            continue;
        end
        colnames{coli} = condition_names{coli-1}{:};
    end
    corresponding_MUs_to_fetch = num2cell(corresponding_MUs_to_fetch);
    corresponding_MUs_to_fetch = cell2table(corresponding_MUs_to_fetch);
    corresponding_MUs_to_fetch.Properties.VariableNames = colnames;
    files_each_condition{nb_of_conditions}.corresponding_MUs_to_fetch = corresponding_MUs_to_fetch;
    files_each_condition{nb_of_conditions}.corresponding_matched_MUs = [];
    files_each_condition{nb_of_conditions}.corresponding_matched_MUs(:,1) = 1:size(corresponding_MUs_to_fetch,1);
    files_each_condition{nb_of_conditions}.corresponding_matched_MUs(:,2) = corresponding_MUs_to_fetch{:,1};
    files_each_condition{nb_of_conditions-1}.corresponding_MUs_to_fetch = corresponding_MUs_to_fetch;
    files_each_condition{nb_of_conditions-1}.corresponding_matched_MUs = [];
    files_each_condition{nb_of_conditions-1}.corresponding_matched_MUs(:,1) = 1:size(corresponding_MUs_to_fetch,1);
    files_each_condition{nb_of_conditions-1}.corresponding_matched_MUs(:,2) = corresponding_MUs_to_fetch{:,1};

    % Add third column to know which ones were included in the FA
    files_each_condition{nb_of_conditions}.corresponding_matched_MUs(:,3) = -2;
    files_each_condition{nb_of_conditions-1}.corresponding_matched_MUs(:,3) = -2;
    %       -2 if never matched
    %       -1 if not matched with plateau
    %       0 if matched with plateau but not included in FA
    %       1 if matched with plateau and included in FA
    for mui = 1:size(corresponding_MUs_to_fetch,1)
        matched_idx_temp = corresponding_MUs_to_fetch{mui,1}; %1st column = matched idx
        if matched_idx_temp > 0
            FA_output_mu_idx_correspondance = find(MUs_clustering_and_correspondance.MU_idx_in_matched_MUs == matched_idx_temp);
            if isempty(FA_output_mu_idx_correspondance)
                % The MU was matched at least once, but not found in the plateau
                % (reference)
                files_each_condition{nb_of_conditions}.corresponding_matched_MUs(mui,3) = -1;
                files_each_condition{nb_of_conditions-1}.corresponding_matched_MUs(mui,3) = -1;
            else
                if MUs_clustering_and_correspondance.Status{FA_output_mu_idx_correspondance} == "matched_discontinuous"
                    % The MU was matched with the plateau, but the MU wasn't
                    % used for the FA
                    files_each_condition{nb_of_conditions}.corresponding_matched_MUs(mui,3) = 0;
                    files_each_condition{nb_of_conditions-1}.corresponding_matched_MUs(mui,3) = 0;
                elseif MUs_clustering_and_correspondance.Status{FA_output_mu_idx_correspondance} == "matched_continuous"
                    % The MU was matched with the plateau, and the MU was used
                    % for the FA
                    files_each_condition{nb_of_conditions}.corresponding_matched_MUs(mui,3) = 1;
                    files_each_condition{nb_of_conditions-1}.corresponding_matched_MUs(mui,3) = 1;
                end
            end
        end
    end
    %

    % Fill the MUs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    time_downsampled_total = 0;
    time_downsampled_total_no_plateau = 0;
    files_each_condition{nb_of_conditions}.time = 0;
    files_each_condition{nb_of_conditions-1}.time = 0;
    for conditioni=1:(nb_of_conditions-2)

        if sum(conditioni == empty_conditions) >= 1
            % Case of an unrepresented condition
            continue;
        end

        time_downsampled_total = time_downsampled_total + size(files_each_condition{conditioni}.smoothed_DR_downsampled,2);
        files_each_condition{nb_of_conditions}.time = files_each_condition{nb_of_conditions}.time +...
            size(files_each_condition{conditioni}.smoothed_DR_concatenated_no_nans,2);
        if conditioni > 1 % for the concatenated file without plateau condition
            time_downsampled_total_no_plateau = time_downsampled_total_no_plateau + ...
                size(files_each_condition{conditioni}.smoothed_DR_downsampled,2);
            files_each_condition{nb_of_conditions-1}.time = files_each_condition{nb_of_conditions-1}.time +...
                size(files_each_condition{conditioni}.smoothed_DR_concatenated_no_nans,2);
        end
    end

    files_each_condition{nb_of_conditions}.time = (1:files_each_condition{nb_of_conditions}.time)./fsamp;
    files_each_condition{nb_of_conditions-1}.time = (1:files_each_condition{nb_of_conditions-1}.time)./fsamp;

    empty_mat_all_matched_MUs_concat = nan(size(corresponding_MUs_to_fetch,1),numel(files_each_condition{nb_of_conditions}.time));
    empty_mat_all_matched_MUs_concat_no_plateau = nan(size(corresponding_MUs_to_fetch,1),...
        numel(files_each_condition{nb_of_conditions-1}.time));
    smoothed_DR_downsampled_all_conditions = nan(size(corresponding_MUs_to_fetch,1),time_downsampled_total);
    smoothed_DR_downsampled_all_conditions_no_plateau = nan(size(corresponding_MUs_to_fetch,1),time_downsampled_total_no_plateau);

    %% Loop through MUs to assign values (all conditions, WITH plateau)

    smoothed_DR_concatenated_no_nans_all_conditions = empty_mat_all_matched_MUs_concat;
    smoothed_DR_concatenated_with_nans_all_conditions = empty_mat_all_matched_MUs_concat;
    files_each_condition{nb_of_conditions}.force = [];
    files_each_condition{nb_of_conditions}.condition_edges = [];
    files_each_condition{nb_of_conditions}.trials_edges = {}; % 1 cell per condition, 1 double per trial
    force_temp = [];

    for mui=1:size(corresponding_MUs_to_fetch,1)
        time_start_condition = 1;
        time_end_condition = 0;
        time_start_condition_downsampled = 1;
        time_end_condition_downsampled = 0;

        for conditioni=1:(nb_of_conditions-2)

            if sum(conditioni == empty_conditions) >= 1
                % Case of an unrepresented condition
                continue;
            end

            time_end_condition = time_start_condition + size(files_each_condition{conditioni}.smoothed_DR_concatenated_no_nans,2) - 1;
            time_end_condition_downsampled = time_start_condition_downsampled + ...
                size(files_each_condition{conditioni}.smoothed_DR_downsampled,2) - 1;

            if mui==1 % Fill window edges (need to be done only for 1st MU) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                files_each_condition{nb_of_conditions}.condition_edges(1:2,conditioni) = [time_start_condition;time_end_condition];
                trials_edges_temp = files_each_condition{conditioni}.window_edges;

                force_temp = [];
                for triali=1:size(trials_edges_temp,1)
                    force_temp = [force_temp,...
                        files_each_condition{conditioni}.force(trials_edges_temp(triali,1):trials_edges_temp(triali,2))];
                end
                trials_edges_temp = trials_edges_temp(:,2) - trials_edges_temp(:,1);
                trials_edges_temp = cumsum(trials_edges_temp)' + (1:numel(trials_edges_temp)); % +1 ; +2 ; +n ... sample for each 1, 2, nth... trial
                trials_edges_temp = trials_edges_temp + (time_start_condition-1);
                files_each_condition{nb_of_conditions}.trials_edges{conditioni} = trials_edges_temp;
                files_each_condition{nb_of_conditions}.force = [files_each_condition{nb_of_conditions}.force,...
                    force_temp];
            end

            %             time_start_condition = time_end_condition + 1;
            %         end

            mu_idx_to_fetch = corresponding_MUs_to_fetch{mui,conditioni+1};

            smoothed_DR_downsampled_all_conditions(mui,time_start_condition_downsampled:time_end_condition_downsampled) = ...
                files_each_condition{conditioni}.smoothed_DR_downsampled(mu_idx_to_fetch,:);

            smoothed_DR_concatenated_no_nans_all_conditions(mui, time_start_condition:time_end_condition) = ...
                files_each_condition{conditioni}.smoothed_DR_concatenated_no_nans(mu_idx_to_fetch,:);
            smoothed_DR_concatenated_with_nans_all_conditions(mui, time_start_condition:time_end_condition) = ...
                files_each_condition{conditioni}.smoothed_DR_concatenated_with_nans(mu_idx_to_fetch,:);

            time_start_condition = time_end_condition + 1;
            time_start_condition_downsampled = time_end_condition_downsampled + 1;

        end % end of "for each condition"
    end % end of "for each MU"

    % Assign values
    files_each_condition{nb_of_conditions}.smoothed_DR_downsampled = smoothed_DR_downsampled_all_conditions;
    files_each_condition{nb_of_conditions}.smoothed_DR_concatenated_no_nans = smoothed_DR_concatenated_no_nans_all_conditions;
    files_each_condition{nb_of_conditions}.smoothed_DR_concatenated_with_nans = smoothed_DR_concatenated_with_nans_all_conditions;

    % Assign window/trials values
    files_each_condition{nb_of_conditions}.window_edges = [];
    files_each_condition{nb_of_conditions}.window_edges_continuous = [];
    files_each_condition{nb_of_conditions}.window_edges_downsampled = [];
    for conditioni=1:(nb_of_conditions-2)
        if sum(conditioni == empty_conditions) >= 1
            % Case of an unrepresented condition
            continue;
        end
        files_each_condition{nb_of_conditions}.window_edges = [...
            files_each_condition{nb_of_conditions}.window_edges;...
            files_each_condition{conditioni}.window_edges];
        temp_window_sample_add = max(files_each_condition{nb_of_conditions}.window_edges_continuous(:));
        if isempty(temp_window_sample_add)
            temp_window_sample_add = 0;
        end
        temp_window_sample_add = temp_window_sample_add + 1;
        files_each_condition{nb_of_conditions}.window_edges_continuous = [...
            files_each_condition{nb_of_conditions}.window_edges_continuous;...
            (files_each_condition{conditioni}.window_edges_continuous) + temp_window_sample_add];
        temp_window_sample_add = max(files_each_condition{nb_of_conditions}.window_edges_downsampled(:));
        if isempty(temp_window_sample_add)
            temp_window_sample_add = 0;
        end
        temp_window_sample_add = temp_window_sample_add;
        files_each_condition{nb_of_conditions}.window_edges_downsampled = [...
            files_each_condition{nb_of_conditions}.window_edges_downsampled;...
            files_each_condition{conditioni}.window_edges_downsampled + temp_window_sample_add];
    end
    % Fill MU data of individual trials
    files_each_condition{nb_of_conditions}.smoothed_DR_per_trial_no_nans = {};
    files_each_condition{nb_of_conditions}.smoothed_DR_per_trial_with_nans = {};
    files_each_condition{nb_of_conditions}.smoothed_DR_per_trial_downsampled = {};
    for triali = 1:size(files_each_condition{nb_of_conditions}.window_edges_continuous,1)
        files_each_condition{nb_of_conditions}.smoothed_DR_per_trial_no_nans{triali} = ...
            files_each_condition{nb_of_conditions}.smoothed_DR_concatenated_no_nans(:,...
            files_each_condition{nb_of_conditions}.window_edges_continuous(triali,1):...
            files_each_condition{nb_of_conditions}.window_edges_continuous(triali,2));
        files_each_condition{nb_of_conditions}.smoothed_DR_per_trial_with_nans{triali} = ...
            files_each_condition{nb_of_conditions}.smoothed_DR_concatenated_with_nans(:,...
            files_each_condition{nb_of_conditions}.window_edges_continuous(triali,1):...
            files_each_condition{nb_of_conditions}.window_edges_continuous(triali,2));
        files_each_condition{nb_of_conditions}.smoothed_DR_per_trial_downsampled{triali} = ...
            files_each_condition{nb_of_conditions}.smoothed_DR_downsampled(:,...
            files_each_condition{nb_of_conditions}.window_edges_downsampled(triali,1):...
            files_each_condition{nb_of_conditions}.window_edges_downsampled(triali,2));
    end

    %% Loop through MUs to assign values (all conditions, WITHOUT plateau)

    smoothed_DR_concatenated_no_nans_all_conditions_no_plateau = empty_mat_all_matched_MUs_concat_no_plateau;
    smoothed_DR_concatenated_with_nans_all_conditions_no_plateau = empty_mat_all_matched_MUs_concat_no_plateau;

    files_each_condition{nb_of_conditions-1}.force = [];
    files_each_condition{nb_of_conditions-1}.condition_edges = [];
    files_each_condition{nb_of_conditions-1}.trials_edges = {}; % 1 cell per condition, 1 double per trial
    force_temp_no_plateau = [];

    for mui=1:size(corresponding_MUs_to_fetch,1)
        time_start_condition = 1;
        time_end_condition = 0;
        time_start_condition_downsampled = 1;
        time_end_condition_downsampled = 0;

        for conditioni=2:(nb_of_conditions-2)

            if sum(conditioni == empty_conditions) >= 1
                % Case of an unrepresented condition
                continue;
            end

            time_end_condition = time_start_condition + size(files_each_condition{conditioni}.smoothed_DR_concatenated_no_nans,2) - 1;
            time_end_condition_downsampled = time_start_condition_downsampled + ...
                size(files_each_condition{conditioni}.smoothed_DR_downsampled,2) - 1;

            if mui==1 % Fill window edges (need to be done only for 1st MU) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                files_each_condition{nb_of_conditions-1}.condition_edges(1:2,conditioni) = [time_start_condition;time_end_condition];
                trials_edges_temp = files_each_condition{conditioni}.window_edges;

                force_temp_no_plateau = [];
                for triali=1:size(trials_edges_temp,1)
                    force_temp_no_plateau = [force_temp_no_plateau,...
                        files_each_condition{conditioni}.force(trials_edges_temp(triali,1):trials_edges_temp(triali,2))];
                end
                trials_edges_temp = trials_edges_temp(:,2) - trials_edges_temp(:,1);
                trials_edges_temp = cumsum(trials_edges_temp)' + (1:numel(trials_edges_temp)); % +1 ; +2 ; +n ... sample for each 1, 2, nth... trial
                trials_edges_temp = trials_edges_temp + (time_start_condition-1);
                files_each_condition{nb_of_conditions-1}.trials_edges{conditioni} = trials_edges_temp;
                files_each_condition{nb_of_conditions-1}.force = [files_each_condition{nb_of_conditions-1}.force,...
                    force_temp_no_plateau];
            end

            mu_idx_to_fetch = corresponding_MUs_to_fetch{mui,conditioni+1};

            smoothed_DR_downsampled_all_conditions_no_plateau(mui,time_start_condition_downsampled:time_end_condition_downsampled) = ...
                files_each_condition{conditioni}.smoothed_DR_downsampled(mu_idx_to_fetch,:);

            smoothed_DR_concatenated_no_nans_all_conditions_no_plateau(mui, time_start_condition:time_end_condition) = ...
                files_each_condition{conditioni}.smoothed_DR_concatenated_no_nans(mu_idx_to_fetch,:);
            smoothed_DR_concatenated_with_nans_all_conditions_no_plateau(mui, time_start_condition:time_end_condition) = ...
                files_each_condition{conditioni}.smoothed_DR_concatenated_with_nans(mu_idx_to_fetch,:);

            time_start_condition = time_end_condition + 1;
            time_start_condition_downsampled = time_end_condition_downsampled + 1;

        end % end of "for each condition"
    end % end of "for each MU"

    % Assign values
    files_each_condition{nb_of_conditions-1}.smoothed_DR_downsampled = smoothed_DR_downsampled_all_conditions_no_plateau;
    files_each_condition{nb_of_conditions-1}.smoothed_DR_concatenated_no_nans = smoothed_DR_concatenated_no_nans_all_conditions_no_plateau;
    files_each_condition{nb_of_conditions-1}.smoothed_DR_concatenated_with_nans = smoothed_DR_concatenated_with_nans_all_conditions_no_plateau;

    % Assign window/trials values
    files_each_condition{nb_of_conditions-1}.window_edges = [];
    files_each_condition{nb_of_conditions-1}.window_edges_continuous = [];
    files_each_condition{nb_of_conditions-1}.window_edges_downsampled = [];
    for conditioni=2:(nb_of_conditions-2)
        if sum(conditioni == empty_conditions) >= 1
            % Case of an unrepresented condition
            continue;
        end
        files_each_condition{nb_of_conditions-1}.window_edges = [...
            files_each_condition{nb_of_conditions-1}.window_edges;...
            files_each_condition{conditioni}.window_edges];
        temp_window_sample_add = max(files_each_condition{nb_of_conditions-1}.window_edges_continuous(:));
        if isempty(temp_window_sample_add)
            temp_window_sample_add = 0;
        end
        temp_window_sample_add = temp_window_sample_add + 1;
        files_each_condition{nb_of_conditions-1}.window_edges_continuous = [...
            files_each_condition{nb_of_conditions-1}.window_edges_continuous;...
            (files_each_condition{conditioni}.window_edges_continuous) + temp_window_sample_add];
        temp_window_sample_add = max(files_each_condition{nb_of_conditions-1}.window_edges_downsampled(:));
        if isempty(temp_window_sample_add)
            temp_window_sample_add = 0;
        end
        temp_window_sample_add = temp_window_sample_add;
        files_each_condition{nb_of_conditions-1}.window_edges_downsampled = [...
            files_each_condition{nb_of_conditions-1}.window_edges_downsampled;...
            files_each_condition{conditioni}.window_edges_downsampled + temp_window_sample_add];
    end
    % Fill MU data of individual trials
    files_each_condition{nb_of_conditions-1}.smoothed_DR_per_trial_no_nans = {};
    files_each_condition{nb_of_conditions-1}.smoothed_DR_per_trial_with_nans = {};
    files_each_condition{nb_of_conditions-1}.smoothed_DR_per_trial_downsampled = {};
    for triali = 1:size(files_each_condition{nb_of_conditions-1}.window_edges_continuous,1)
        files_each_condition{nb_of_conditions-1}.smoothed_DR_per_trial_no_nans{triali} = ...
            files_each_condition{nb_of_conditions-1}.smoothed_DR_concatenated_no_nans(:,...
            files_each_condition{nb_of_conditions-1}.window_edges_continuous(triali,1):...
            files_each_condition{nb_of_conditions-1}.window_edges_continuous(triali,2));
        files_each_condition{nb_of_conditions-1}.smoothed_DR_per_trial_with_nans{triali} = ...
            files_each_condition{nb_of_conditions-1}.smoothed_DR_concatenated_with_nans(:,...
            files_each_condition{nb_of_conditions-1}.window_edges_continuous(triali,1):...
            files_each_condition{nb_of_conditions-1}.window_edges_continuous(triali,2));
        files_each_condition{nb_of_conditions-1}.smoothed_DR_per_trial_downsampled{triali} = ...
            files_each_condition{nb_of_conditions-1}.smoothed_DR_downsampled(:,...
            files_each_condition{nb_of_conditions-1}.window_edges_downsampled(triali,1):...
            files_each_condition{nb_of_conditions-1}.window_edges_downsampled(triali,2));
    end

    %% Display full concatenated signal, to check if everything looks right
    close(figure(999))
    figure(999)
    colors_concat_temp = turbo(size(corresponding_MUs_to_fetch,1));
    mean_concat_DR = [];
    for mui=1:size(corresponding_MUs_to_fetch,1)
        mean_concat_DR(mui) = mean(files_each_condition{nb_of_conditions}.smoothed_DR_concatenated_no_nans(mui,:),'omitnan');
    end
    plot(files_each_condition{nb_of_conditions}.time,...
        files_each_condition{nb_of_conditions}.force ./ max(files_each_condition{nb_of_conditions}.force) .* max(mean_concat_DR),...
        'color',[.25,.25,.25],'LineWidth',2);
    hold on
    [~,order_by_mean_DR] = sort(mean_concat_DR,'ascend');
    for mui=1:size(corresponding_MUs_to_fetch,1)
        mu_to_plot_idx = order_by_mean_DR(mui);
        plot(files_each_condition{nb_of_conditions}.time,...
            files_each_condition{nb_of_conditions}.smoothed_DR_concatenated_no_nans(mu_to_plot_idx,:), ...
            'color',colors_concat_temp(mui,:),'linew',1);
    end
    for conditioni=1:(nb_of_conditions-2)
        if sum(conditioni == empty_conditions) >= 1
            % Case of an unrepresented condition
            continue;
        end
        for triali = 1:numel(files_each_condition{nb_of_conditions}.trials_edges{conditioni})
            line([files_each_condition{nb_of_conditions}.trials_edges{conditioni}(triali),...
                files_each_condition{nb_of_conditions}.trials_edges{conditioni}(triali)] ./ fsamp,...
                ylim(),'Color','k','LineStyle','--','LineWidth',1.5);
        end
        line([files_each_condition{nb_of_conditions}.condition_edges(2,conditioni),...
            files_each_condition{nb_of_conditions}.condition_edges(2,conditioni)] ./ fsamp,...
            ylim(),'Color','k','LineWidth',2);
    end
    title("MUs matched in all files, conctenated for all conditions");
    xlabel("Time (s)");
    ylabel("Discharge rate (spike/s) | Force (no unit ; arbitrarily scaled for visibility)");

    saveas(gcf(),strcat(output_folder,'/matched_MUs_all_conditions_concat.png'));

    % SAME, WITHOUT PLATEAU
    close(figure(9999))
    figure(9999)
    colors_concat_temp = turbo(size(corresponding_MUs_to_fetch,1));
    mean_concat_DR = [];
    for mui=1:size(corresponding_MUs_to_fetch,1)
        mean_concat_DR(mui) = mean(files_each_condition{nb_of_conditions-1}.smoothed_DR_concatenated_no_nans(mui,:),'omitnan');
    end
    plot(files_each_condition{nb_of_conditions-1}.time,...
        files_each_condition{nb_of_conditions-1}.force ./ max(files_each_condition{nb_of_conditions-1}.force) .* max(mean_concat_DR),...
        'color',[.25,.25,.25],'LineWidth',2);
    hold on
    [~,order_by_mean_DR] = sort(mean_concat_DR,'ascend');
    for mui=1:size(corresponding_MUs_to_fetch,1)
        mu_to_plot_idx = order_by_mean_DR(mui);
        plot(files_each_condition{nb_of_conditions-1}.time,...
            files_each_condition{nb_of_conditions-1}.smoothed_DR_concatenated_no_nans(mu_to_plot_idx,:), ...
            'color',colors_concat_temp(mui,:),'linew',1);
    end
    for conditioni=1:(nb_of_conditions-2)
        if sum(conditioni == empty_conditions) >= 1
            % Case of an unrepresented condition
            continue;
        end
        if isempty(files_each_condition{nb_of_conditions-1}.trials_edges{conditioni})
            continue;
        end
        for triali = 1:numel(files_each_condition{nb_of_conditions-1}.trials_edges{conditioni})
            line([files_each_condition{nb_of_conditions-1}.trials_edges{conditioni}(triali),...
                files_each_condition{nb_of_conditions-1}.trials_edges{conditioni}(triali)] ./ fsamp,...
                ylim(),'Color','k','LineStyle','--','LineWidth',1.5);
        end
        line([files_each_condition{nb_of_conditions-1}.condition_edges(2,conditioni),...
            files_each_condition{nb_of_conditions-1}.condition_edges(2,conditioni)] ./ fsamp,...
            ylim(),'Color','k','LineWidth',2);
    end
    title("MUs matched in all files, conctenated for all conditions (without reference condition)");
    xlabel("Time (s)");
    ylabel("Discharge rate (spike/s) | Force (no unit ; arbitrarily scaled for visibility)");

    saveas(gcf(),strcat(output_folder,'/matched_MUs_all_conditions_concat_no_plateau.png'));

end % end of "if add a condition with all conditions concatenated"

%% COMPUTE FLEXIBILITY OF ENTIRE POPULATION (as displacement and dispersion)

% Careful, the "reponsible units" values are not the matched MU indexes,
% but the indexes in the file

if compute_flexibility_for_entire_population

    max_flexibility_overall = 0; % this is used only to plot dispersions on the same scale later in the script

    for conditioni=1:nb_of_conditions %=2

        %% Initiliaze values

        files_each_condition{conditioni}.flexibility = struct();
        files_each_condition{conditioni}.flexibility.nb_of_MUs_in_pop_for_flexibility = 0;
        files_each_condition{conditioni}.flexibility.nb_of_MUs_in_pop_for_flexibility_per_trial = {};

        % dispersion
        files_each_condition{conditioni}.flexibility.population_dispersions = nan;
        files_each_condition{conditioni}.flexibility.population_dispersions_per_trial = {};
        % dispersion - Max
        files_each_condition{conditioni}.flexibility.max_population_dispersion = nan;
        files_each_condition{conditioni}.flexibility.max_population_dispersion_normalized_by_MU_nb = nan;
        files_each_condition{conditioni}.flexibility.max_population_dispersion_normalized_by_mean_DR = nan;
        files_each_condition{conditioni}.flexibility.max_population_dispersion_per_trial = {};
        files_each_condition{conditioni}.flexibility.max_population_dispersion_normalized_by_MU_nb_per_trial = {};
        files_each_condition{conditioni}.flexibility.max_population_dispersion_normalized_by_mean_DR_per_trial = {};
        files_each_condition{conditioni}.flexibility.max_population_dispersion_mean_per_trial = nan;
        files_each_condition{conditioni}.flexibility.max_population_dispersion_normalized_by_MU_nb_mean_per_trial = nan;
        files_each_condition{conditioni}.flexibility.max_population_dispersion_normalized_by_mean_DR_mean_per_trial = nan;
        % dispersion - Mean
        files_each_condition{conditioni}.flexibility.mean_population_dispersion = nan;
        files_each_condition{conditioni}.flexibility.mean_population_dispersion_normalized_by_MU_nb = nan;
        files_each_condition{conditioni}.flexibility.mean_population_dispersion_normalized_by_mean_DR = nan;
        files_each_condition{conditioni}.flexibility.mean_population_dispersion_per_trial = {};
        files_each_condition{conditioni}.flexibility.mean_population_dispersion_normalized_by_MU_nb_per_trial = {};
        files_each_condition{conditioni}.flexibility.mean_population_dispersion_normalized_by_mean_DR_per_trial = {};
        files_each_condition{conditioni}.flexibility.mean_population_dispersion_mean_per_trial = nan;
        files_each_condition{conditioni}.flexibility.mean_population_dispersion_normalized_by_MU_nb_mean_per_trial = nan;
        files_each_condition{conditioni}.flexibility.mean_population_dispersion_normalized_by_mean_DR_mean_per_trial = nan;

        % displacement
        files_each_condition{conditioni}.flexibility.population_displacements = nan;
        files_each_condition{conditioni}.flexibility.population_displacements_per_trial = {};
        files_each_condition{conditioni}.flexibility.population_displacements_responsible_units = nan;
        files_each_condition{conditioni}.flexibility.population_displacements_responsible_units_per_trial = {};
        % displacement - Max
        files_each_condition{conditioni}.flexibility.max_population_displacement = nan;
        files_each_condition{conditioni}.flexibility.max_population_displacement_normalized_by_MU_nb = nan;
        files_each_condition{conditioni}.flexibility.max_population_displacement_normalized_by_mean_DR = nan;
        files_each_condition{conditioni}.flexibility.max_population_displacement_per_trial = {};
        files_each_condition{conditioni}.flexibility.max_population_displacement_normalized_by_MU_nb_per_trial = {};
        files_each_condition{conditioni}.flexibility.max_population_displacement_normalized_by_mean_DR_per_trial = {};
        files_each_condition{conditioni}.flexibility.max_population_displacement_mean_per_trial = nan;
        files_each_condition{conditioni}.flexibility.max_population_displacement_normalized_by_MU_nb_mean_per_trial = nan;
        files_each_condition{conditioni}.flexibility.max_population_displacement_normalized_by_mean_DR_mean_per_tria = nan; %not 'trial' because of character limit
        % displacement - Mean
        files_each_condition{conditioni}.flexibility.mean_population_displacement = nan;
        files_each_condition{conditioni}.flexibility.mean_population_displacement_normalized_by_MU_nb = nan;
        files_each_condition{conditioni}.flexibility.mean_population_displacement_normalized_by_mean_DR = nan;
        files_each_condition{conditioni}.flexibility.mean_population_displacement_per_trial = {};
        files_each_condition{conditioni}.flexibility.mean_population_displacement_normalized_by_MU_nb_per_trial = {};
        files_each_condition{conditioni}.flexibility.mean_population_displacement_normalized_by_mean_DR_per_trial = {};
        files_each_condition{conditioni}.flexibility.mean_population_displacement_mean_per_trial = nan;
        files_each_condition{conditioni}.flexibility.mean_population_displacement_normalized_by_MU_nb_mean_per_trial = nan;
        files_each_condition{conditioni}.flexibility.mean_population_displacement_normalized_by_mean_DR_mean_per_tri = nan; %not 'trial' because of character limit


        %% end of initializing values

        if sum(conditioni == empty_conditions) >= 1
            % DO NOTHING = condition not represented in this subject
            save(strcat(output_folder,"/condition_",num2str(conditioni),"_is_empty.txt"),"empty_conditions");

        else
            % Normal case of a represented condition

            if downsampling
                DRs_to_use = files_each_condition{conditioni}.smoothed_DR_downsampled;
                DRs_to_use_per_trial = files_each_condition{conditioni}.smoothed_DR_per_trial_downsampled;
                fsamp_to_use = fsamp/downsampling_factor;
            else
                DRs_to_use = files_each_condition{conditioni}.smoothed_DR_concatenated_no_nans;
                DRs_to_use_per_trial = files_each_condition{conditioni}.smoothed_DR_per_trial_no_nans;
                fsamp_to_use = fsamp;
            end

            %% FLEXIBILITY PER TRIAL
            for triali = 1:numel(DRs_to_use_per_trial)
                MUs_DRs_as_psth = {};
                MUs_DRs_as_psth{1} = DRs_to_use_per_trial{triali};

                files_each_condition{conditioni}.flexibility.nb_of_MUs_in_pop_for_flexibility_per_trial{triali} = ...
                    size(DRs_to_use_per_trial{triali},1);

                % DISPLACEMENT
                displacements_temp = [];
                responsibleUnits = [];
                [displacements_temp, ~, responsibleUnits] = ...
                    mudisplacement(MUs_DRs_as_psth,fsamp_to_use);

                files_each_condition{conditioni}.flexibility.population_displacements_per_trial{triali} = ...
                    displacements_temp;
                files_each_condition{conditioni}.flexibility.population_displacements_responsible_units_per_trial{triali} = ...
                    responsibleUnits;
                % Displacement MAX
                files_each_condition{conditioni}.flexibility.max_population_displacement_per_trial{triali} = ...
                    max(displacements_temp);
                files_each_condition{conditioni}.flexibility.max_population_displacement_normalized_by_MU_nb_per_trial{triali} = ...
                    max(displacements_temp)./size(MUs_DRs_as_psth{1},1);
                files_each_condition{conditioni}.flexibility.max_population_displacement_normalized_by_mean_DR_per_trial{triali} = ...
                    max(displacements_temp)./mean(mean(MUs_DRs_as_psth{1},2,'omitnan'),'omitnan');
                % Displacement MEAN
                files_each_condition{conditioni}.flexibility.mean_population_displacement_per_trial{triali} = ...
                    mean(displacements_temp,'omitnan');
                files_each_condition{conditioni}.flexibility.mean_population_displacement_normalized_by_MU_nb_per_trial{triali} = ...
                    mean(displacements_temp,'omitnan')./size(MUs_DRs_as_psth{1},1);
                files_each_condition{conditioni}.flexibility.mean_population_displacement_normalized_by_mean_DR_per_trial{triali} = ...
                    mean(displacements_temp,'omitnan')./mean(mean(MUs_DRs_as_psth{1},2,'omitnan'),'omitnan');

                % DISPERSION
                population_dispersion = [];
                [population_dispersion, ~, ~] = ...
                    mnpdispersion(MUs_DRs_as_psth,fsamp_to_use,'lagType','global'); % 'local' is faster but can cause inconsistencies

                files_each_condition{conditioni}.flexibility.population_dispersions_per_trial{triali} = population_dispersion;
                % Dispersion MAX
                files_each_condition{conditioni}.flexibility.max_population_dispersion_per_trial{triali} = ...
                    max(population_dispersion);
                files_each_condition{conditioni}.flexibility.max_population_dispersion_normalized_by_MU_nb_per_trial{triali} = ...
                    max(population_dispersion)./size(MUs_DRs_as_psth{1},1);
                files_each_condition{conditioni}.flexibility.max_population_dispersion_normalized_by_mean_DR_per_trial{triali} = ...
                    max(population_dispersion)./mean(mean(MUs_DRs_as_psth{1},2,'omitnan'),'omitnan');
                % Dispersion MEAN
                files_each_condition{conditioni}.flexibility.mean_population_dispersion_per_trial{triali} = ...
                    mean(population_dispersion,'omitnan');
                files_each_condition{conditioni}.flexibility.mean_population_dispersion_normalized_by_MU_nb_per_trial{triali} = ...
                    mean(population_dispersion,'omitnan')./size(MUs_DRs_as_psth{1},1);
                files_each_condition{conditioni}.flexibility.mean_population_dispersion_normalized_by_mean_DR_per_trial{triali} = ...
                    mean(population_dispersion,'omitnan')./mean(mean(MUs_DRs_as_psth{1},2,'omitnan'),'omitnan');

            end % end of "for each trial"

            % MEAN PER TRIAL VALUES
            % Displacement - MAX
            files_each_condition{conditioni}.flexibility.max_population_displacement_mean_per_trial = ...
                mean([files_each_condition{conditioni}.flexibility.max_population_displacement_per_trial{:}],'omitnan');
            files_each_condition{conditioni}.flexibility.max_population_displacement_normalized_by_MU_nb_mean_per_trial = ...
                mean([files_each_condition{conditioni}.flexibility.max_population_displacement_normalized_by_MU_nb_per_trial{:}],'omitnan');
            files_each_condition{conditioni}.flexibility.max_population_displacement_normalized_by_mean_DR_mean_per_tria = ... %not 'trial' because of character limit
                mean([files_each_condition{conditioni}.flexibility.max_population_displacement_normalized_by_mean_DR_per_trial{:}],'omitnan');
            % Displacement - MEAN
            files_each_condition{conditioni}.flexibility.mean_population_displacement_mean_per_trial = ...
                mean([files_each_condition{conditioni}.flexibility.mean_population_displacement_per_trial{:}],'omitnan');
            files_each_condition{conditioni}.flexibility.mean_population_displacement_normalized_by_MU_nb_mean_per_trial = ...
                mean([files_each_condition{conditioni}.flexibility.mean_population_displacement_normalized_by_MU_nb_per_trial{:}],'omitnan');
            files_each_condition{conditioni}.flexibility.mean_population_displacement_normalized_by_mean_DR_mean_per_tri = ... %not 'trial' because of character limit
                mean([files_each_condition{conditioni}.flexibility.mean_population_displacement_normalized_by_mean_DR_per_trial{:}],'omitnan');

            % Dispersion - MAX
            files_each_condition{conditioni}.flexibility.max_population_dispersion_mean_per_trial = ...
                mean([files_each_condition{conditioni}.flexibility.max_population_dispersion_per_trial{:}],'omitnan');
            files_each_condition{conditioni}.flexibility.max_population_dispersion_normalized_by_MU_nb_mean_per_trial = ...
                mean([files_each_condition{conditioni}.flexibility.max_population_dispersion_normalized_by_MU_nb_per_trial{:}],'omitnan');
            files_each_condition{conditioni}.flexibility.max_population_dispersion_normalized_by_mean_DR_mean_per_trial = ...
                mean([files_each_condition{conditioni}.flexibility.max_population_dispersion_normalized_by_mean_DR_per_trial{:}],'omitnan');
            % Dispersion - MEAN
            files_each_condition{conditioni}.flexibility.mean_population_dispersion_mean_per_trial = ...
                mean([files_each_condition{conditioni}.flexibility.mean_population_dispersion_per_trial{:}],'omitnan');
            files_each_condition{conditioni}.flexibility.mean_population_dispersion_normalized_by_MU_nb_mean_per_trial = ...
                mean([files_each_condition{conditioni}.flexibility.mean_population_dispersion_normalized_by_MU_nb_per_trial{:}],'omitnan');
            files_each_condition{conditioni}.flexibility.mean_population_dispersion_normalized_by_mean_DR_mean_per_trial = ... %not 'trial' because of character limit
                mean([files_each_condition{conditioni}.flexibility.mean_population_dispersion_normalized_by_mean_DR_per_trial{:}],'omitnan');

            %% FLEXIBILITY FOR TRIALS CONCATENATED

            MUs_DRs_as_psth = {};
            MUs_DRs_as_psth{1} = DRs_to_use;

            files_each_condition{conditioni}.flexibility.nb_of_MUs_in_pop_for_flexibility = ...
                size(MUs_DRs_as_psth{1},1);

            % DISPLACEMENT
            displacements_temp = [];
            responsibleUnits = [];
            [displacements_temp, tau, responsibleUnits] = ...
                mudisplacement(MUs_DRs_as_psth,fsamp_to_use);

            files_each_condition{conditioni}.flexibility.population_displacements = displacements_temp;
            files_each_condition{conditioni}.flexibility.population_displacements_responsible_units = responsibleUnits;
            % Displacement max
            files_each_condition{conditioni}.flexibility.max_population_displacement = ...
                max(displacements_temp);
            files_each_condition{conditioni}.flexibility.max_population_displacement_normalized_by_MU_nb = ...
                max(displacements_temp)./size(MUs_DRs_as_psth{1},1);
            files_each_condition{conditioni}.flexibility.max_population_displacement_normalized_by_mean_DR = ...
                max(displacements_temp)./mean(mean(MUs_DRs_as_psth{1},2,'omitnan'),'omitnan');
            % Displacement mean
            files_each_condition{conditioni}.flexibility.mean_population_displacement = ...
                mean(displacements_temp,'omitnan');
            files_each_condition{conditioni}.flexibility.mean_population_displacement_normalized_by_MU_nb = ...
                mean(displacements_temp,'omitnan')./size(MUs_DRs_as_psth{1},1);
            files_each_condition{conditioni}.flexibility.mean_population_displacement_normalized_by_mean_DR = ...
                mean(displacements_temp,'omitnan')./mean(mean(MUs_DRs_as_psth{1},2,'omitnan'),'omitnan');

            % DISPERSION
            population_dispersion = [];
            [population_dispersion, ~, ~] = ...
                mnpdispersion(MUs_DRs_as_psth,fsamp_to_use,'lagType','global'); % 'local' is faster but can cause inconsistencies

            files_each_condition{conditioni}.flexibility.population_dispersions = population_dispersion;
            % Dispersion max
            files_each_condition{conditioni}.flexibility.max_population_dispersion = ...
                max(population_dispersion);
            files_each_condition{conditioni}.flexibility.max_population_dispersion_normalized_by_MU_nb = ...
                max(population_dispersion)./size(MUs_DRs_as_psth{1},1);
            files_each_condition{conditioni}.flexibility.max_population_dispersion_normalized_by_mean_DR = ...
                max(population_dispersion)./mean(mean(MUs_DRs_as_psth{1},2,'omitnan'),'omitnan');
            % Dispersion mean
            files_each_condition{conditioni}.flexibility.mean_population_dispersion = ...
                mean(population_dispersion,'omitnan');
            files_each_condition{conditioni}.flexibility.mean_population_dispersion_normalized_by_MU_nb = ...
                mean(population_dispersion,'omitnan')./size(MUs_DRs_as_psth{1},1);
            files_each_condition{conditioni}.flexibility.mean_population_dispersion_normalized_by_mean_DR = ...
                mean(population_dispersion,'omitnan')./mean(mean(MUs_DRs_as_psth{1},2,'omitnan'),'omitnan');


        end % End of "if Case of an unrepresented VS represented condition"

        vars_to_save = {};
        % Get field names
        fieldNames = fieldnames(files_each_condition{conditioni}.flexibility);
        % Assign the field names to workspace variables
        for fieldi=1:numel(fieldNames)
            assignin('base', fieldNames{fieldi}, ...
                files_each_condition{conditioni}.flexibility.(fieldNames{fieldi})); % assigning variable to the worksapce
        end
        vars_to_save = fieldNames;
        save(strcat(output_folder,"/",condition_names{conditioni},"/displacement_dispersion_output_trials_concatenated.mat"),...
            vars_to_save{:});

    end % end of "for each condition"

end %end of "if compute flexibility for enire population

%%
if ~compute_flexibility_per_MU_pair
    % If the condition is met, stop executing the script
    return;
end

%%
% COMPUTE FLEXIBILTY PER MU PAIR

for conditioni=1:nb_of_conditions %1:nb_of_conditions

    progressBarCondition = waitbar(conditioni/nb_of_conditions, ...
        strcat("Flexibility for each MU pair ; condition ",num2str(conditioni),"/", ...
        num2str(nb_of_conditions)));

    % Create empty table of flexibility per MU pair for each condition
    colNames = {'MU_x','MU_x_matched','MU_x_FA_status','MU_y','MU_y_matched','MU_y_FA_status',...
        'MU_x_muscle','MU_y_muscle','MU_pair_muscles', ...
        'max_dispersion','mean_dispersion','mean_dispersion_normalized',...
        'max_dispersion_trials_mean','mean_dispersion_trials_mean','mean_dispersion_normalized_trials_mean',...
        'max_displacement','mean_displacement','max_displacement_trials_mean','mean_displacement_trials_mean',...
        'cross_corr','cross_corr_lag','cross_corr_trials_mean','cross_corr_lag_trials_mean'};
    nb_of_cols = numel(colNames);
    empty_row = {};
    for coli=1:nb_of_cols
        empty_row{coli} = {nan};
    end
    flexibility_between_MU_pairs = table(empty_row{:});
    flexibility_between_MU_pairs.Properties.VariableNames = colNames;

    if sum(conditioni == empty_conditions) >= 1
        % Case of an unrepresented condition
        files_each_condition{conditioni}.flexibility_MU_pairs = flexibility_between_MU_pairs;
    else
        % Normal case of a represented condition


        DRs_to_use = files_each_condition{conditioni}.smoothed_DR_concatenated_no_nans;
        DRs_to_use_per_trial = {};
        for triali=1:numel(files_each_condition{conditioni}.smoothed_DR_per_trial_no_nans)
            DRs_to_use_per_trial{triali} = files_each_condition{conditioni}.smoothed_DR_per_trial_no_nans{triali};
        end
        if downsampling && perform_displacement_calculation_per_mu_pair
            displacement_per_mu_pair = true;
            DR_downsampled_to_use = files_each_condition{conditioni}.smoothed_DR_downsampled;
            DR_downsampled_to_use_per_trial = {};
            for triali=1:numel(files_each_condition{conditioni}.smoothed_DR_per_trial_downsampled)
                DR_downsampled_to_use_per_trial{triali} = files_each_condition{conditioni}.smoothed_DR_per_trial_downsampled{triali};
            end
            fsamp_to_use = fsamp/downsampling_factor;
        else
            displacement_per_mu_pair = false;
        end
        pair_iter = 0;
        nb_MUs = size(DRs_to_use,1); %5 %for testing purpose

        MUs_correspondance = files_each_condition{conditioni}.corresponding_matched_MUs;
        % Assign muscle names to each MU
        MUs_muscles = {};
        for mui = 1:size(MUs_correspondance,1)
            % if "all conditions concatenated", we do not have the grid
            if (all_conditions_concatenated) && (conditioni > (nb_of_conditions - 2))
               % find the corresponding grid in MUs matched
               MUs_muscles{mui} = convertCharsToStrings(MUs_matched{MUs_correspondance(mui,2),"Muscle"}{:});
            else
                grid_temp = floor(MUs_correspondance(mui,1)./100);
                MUs_muscles{mui} = muscles_per_grid{grid_temp};
            end

        end

        if display_figures_dispersion && save_figures_dispersion
            new_folder_trials_concat = strcat(output_folder,"/dispersion_figures/",condition_names{conditioni});
            mkdir(new_folder_trials_concat);
            new_folder_single_trial = {};
            for triali=1:numel(DRs_to_use_per_trial)
                new_folder_single_trial{triali} = strcat(output_folder,"/dispersion_figures/",condition_names{conditioni},"/trial_",num2str(triali));
                mkdir(new_folder_single_trial{triali});
            end
        end

        nb_iterations = nchoosek(nb_MUs,2);
        progressBarMUPair = waitbar(0, ...
            strcat("MU pair ",num2str(pair_iter),"/", ...
            num2str(nb_iterations)));

        for mu_xi=1:(nb_MUs-1) %3 %(nb_MUs-1)
            for mu_yi=mu_xi+1:nb_MUs %4 %nb_MUs
                pair_iter = pair_iter+1;

                waitbar(pair_iter/nb_iterations, progressBarMUPair, ...
                    strcat("MU pair ",num2str(pair_iter),"/", ...
                    num2str(nb_iterations)));

                if display_figures_dispersion
                    close(figure(conditioni*1e6 + 1e4*mu_xi + 1e2*mu_yi));
                    for triali=1:numel(DRs_to_use_per_trial)
                        close(figure(conditioni*1e6 + 1e4*mu_xi + 1e2*mu_yi + triali));
                    end
                end

                %% DISPERSION
                % Actual MU pair dispersion calculation
                % Trials concatenated
                [~ , ~ , ~, dispersion_temp , ~] = Dispersion_for_MU_Pairs_with_timelag(  ...
                    DRs_to_use(mu_xi,:), ...
                    DRs_to_use(mu_yi,:), ...
                    fsamp, 1, 0, 0, display_figures_dispersion, conditioni*1e6 + 1e4*mu_xi + 1e2*mu_yi);
                % Each trial separately
                dispersion_temp_single_trial = {};
                for triali=1:numel(DRs_to_use_per_trial)
                    [~ , ~ , ~, dispersion_temp_single_trial{triali} , ~] = Dispersion_for_MU_Pairs_with_timelag(  ...
                        DRs_to_use_per_trial{triali}(mu_xi,:), ...
                        DRs_to_use_per_trial{triali}(mu_yi,:), ...
                        fsamp, 1, 0, 0, display_figures_dispersion, conditioni*1e6 + 1e4*mu_xi + 1e2*mu_yi + triali);
                end

                if display_figures_dispersion && save_figures_dispersion
                    saveas(figure(conditioni*1e6 + 1e4*mu_xi + 1e2*mu_yi), ...
                        strcat(new_folder_trials_concat,"/MU pair (trials concatenated) - MU#", ...
                        num2str(MUs_correspondance(mu_xi,2)), " (in file MU#", num2str(mu_xi), ...
                        ") vs MU#", ...
                        num2str(MUs_correspondance(mu_yi,2)), " (in file MU#", num2str(mu_yi),...
                        ").png")...
                        );
                    for triali=1:numel(DRs_to_use_per_trial)
                        saveas(figure(conditioni*1e6 + 1e4*mu_xi + 1e2*mu_yi + triali), ...
                            strcat(new_folder_single_trial{triali},"/MU pair (trial #",num2str(triali),...
                            ") - MU#", ...
                            num2str(MUs_correspondance(mu_xi,2)), " (in file MU#", num2str(mu_xi), ...
                            ") vs MU#", ...
                            num2str(MUs_correspondance(mu_yi,2)), " (in file MU#", num2str(mu_yi),...
                            ").png")...
                            );
                    end
                end
                flexibility_between_MU_pairs{pair_iter,:} = empty_row{:};
                flexibility_between_MU_pairs{pair_iter,'MU_x'} = {mu_xi};
                flexibility_between_MU_pairs{pair_iter,'MU_x_matched'} = {...
                    MUs_correspondance(mu_xi,2)};
                flexibility_between_MU_pairs{pair_iter,'MU_x_FA_status'} = {...
                    MUs_correspondance(mu_xi,3)};
                flexibility_between_MU_pairs{pair_iter,'MU_y'} = {mu_yi};
                flexibility_between_MU_pairs{pair_iter,'MU_y_matched'} = {...
                    MUs_correspondance(mu_yi,2)};
                flexibility_between_MU_pairs{pair_iter,'MU_y_FA_status'} = {...
                    MUs_correspondance(mu_yi,3)};

                flexibility_between_MU_pairs{pair_iter,'MU_x_muscle'} = {...
                    MUs_muscles{mu_xi}};
                flexibility_between_MU_pairs{pair_iter,'MU_y_muscle'} = {...
                    MUs_muscles{mu_yi}};
                % put into alphabetical order (so that VL-RF for example
                % is the same as RF-VL)
                muscles_alphabet_order = sort({MUs_muscles{mu_xi}{:},MUs_muscles{mu_yi}{:}});
                flexibility_between_MU_pairs{pair_iter,'MU_pair_muscles'} = {...
                    strcat(muscles_alphabet_order{1},"-",muscles_alphabet_order{2})};

                flexibility_between_MU_pairs{pair_iter,'max_dispersion'} = {max([dispersion_temp{:}])};
                flexibility_between_MU_pairs{pair_iter,'mean_dispersion'} = {mean([dispersion_temp{:}],'omitnan')};
                flexibility_between_MU_pairs{pair_iter,'mean_dispersion_normalized'} = ...
                    {mean([dispersion_temp{:}],'omitnan')/numel(dispersion_temp)};
                % per trial
                max_dispersion_temp_per_trial = [];
                mean_dispersion_temp_per_trial = [];
                mean_dispersion_normalized_temp_per_trial = [];
                for triali=1:numel(DRs_to_use_per_trial)
                    if ~isempty(dispersion_temp_single_trial{triali})
                        max_dispersion_temp_per_trial(triali) = ...
                            max([dispersion_temp_single_trial{triali}{:}]);
                        mean_dispersion_temp_per_trial(triali) = ...
                            mean([dispersion_temp_single_trial{triali}{:}]);
                        mean_dispersion_normalized_temp_per_trial(triali) = ...
                            mean([dispersion_temp_single_trial{triali}{:}])/numel(dispersion_temp_single_trial{triali});
                    else
                        max_dispersion_temp_per_trial(triali) = nan;
                        mean_dispersion_temp_per_trial(triali) = nan;
                        mean_dispersion_normalized_temp_per_trial(triali) = nan;
                    end
                end
                flexibility_between_MU_pairs{pair_iter,'max_dispersion_trials_mean'} = {mean(max_dispersion_temp_per_trial,'omitnan')};
                flexibility_between_MU_pairs{pair_iter,'mean_dispersion_trials_mean'} = {mean(mean_dispersion_temp_per_trial,'omitnan')};
                flexibility_between_MU_pairs{pair_iter,'mean_dispersion_normalized_trials_mean'} = {mean(mean_dispersion_normalized_temp_per_trial,'omitnan')};

                %% DISPLACEMENT
                if displacement_per_mu_pair
                    MU_pair_DRs_as_psth = {};
                    MU_pair_DRs_as_psth{1} = [DR_downsampled_to_use(mu_xi,:); ...
                        DR_downsampled_to_use(mu_yi,:)];

                    [displacement_mupair_temp, ~, ~] = mudisplacement(MU_pair_DRs_as_psth,fsamp_to_use);
                    flexibility_between_MU_pairs{pair_iter,'max_displacement'} = {max(displacement_mupair_temp)};
                    flexibility_between_MU_pairs{pair_iter,'mean_displacement'} = {mean(displacement_mupair_temp)};

                    % per trial
                    displacements_temp_single_trial = {};
                    for triali=1:numel(DR_downsampled_to_use_per_trial)
                        MU_pair_DRs_as_psth = {};
                        MU_pair_DRs_as_psth{1} = [DR_downsampled_to_use_per_trial{triali}(mu_xi,:); ...
                            DR_downsampled_to_use_per_trial{triali}(mu_yi,:)];
                        [displacements_temp_single_trial{triali}, ~, ~] = mudisplacement(MU_pair_DRs_as_psth,fsamp_to_use);
                    end

                    max_displacements_temp_per_trial = [];
                    mean_displacements_temp_per_trial = [];
                    for triali=1:numel(DR_downsampled_to_use_per_trial)
                        if ~isempty(displacements_temp_single_trial{triali})
                            max_displacements_temp_per_trial(triali) = ...
                                max([displacements_temp_single_trial{triali}]);
                            mean_displacements_temp_per_trial(triali) = ...
                                mean([displacements_temp_single_trial{triali}]);
                        else
                            max_displacements_temp_per_trial(triali) = nan;
                            mean_displacements_temp_per_trial(triali) = nan;
                        end
                    end
                    flexibility_between_MU_pairs{pair_iter,'max_displacement_trials_mean'} = {mean(max_displacements_temp_per_trial,'omitnan')};
                    flexibility_between_MU_pairs{pair_iter,'mean_displacement_trials_mean'} = {mean(mean_displacements_temp_per_trial,'omitnan')};


                else
                    flexibility_between_MU_pairs{pair_iter,'max_displacement'} = {nan};
                    flexibility_between_MU_pairs{pair_iter,'mean_displacement'} = {nan};
                end

                %% CROSS-CORREL
                muX_for_cross_correl = DRs_to_use(mu_xi,:)';
                muX_for_cross_correl = muX_for_cross_correl - mean(muX_for_cross_correl,'omitnan');
                muY_for_cross_correl = DRs_to_use(mu_yi,:)';
                muY_for_cross_correl = muY_for_cross_correl - mean(muY_for_cross_correl,'omitnan');
                [r_temp, lags] = xcorr(muX_for_cross_correl,muY_for_cross_correl,...
                    round(fsamp*0.1),"normalized");
                [max_corr_value, max_corr_lag] = max(r_temp);
                % Lag of first MU relative to the other (positive values =
                % make the MU DR happen later)
                max_corr_lag = (max_corr_lag - round(length(lags)/2) ) /fsamp ;

                flexibility_between_MU_pairs{pair_iter,'cross_corr'} = {max_corr_value};
                flexibility_between_MU_pairs{pair_iter,'cross_corr_lag'} = {max_corr_lag};

                % per trial
                correl_values_temp_per_trial = [];
                correl_lags_temp_per_trial = [];
                for triali=1:numel(DRs_to_use_per_trial)
                    muX_for_cross_correl = DRs_to_use_per_trial{triali}(mu_xi,:)';
                    muX_for_cross_correl = muX_for_cross_correl - mean(muX_for_cross_correl,'omitnan');
                    muY_for_cross_correl = DRs_to_use_per_trial{triali}(mu_yi,:)';
                    muY_for_cross_correl = muY_for_cross_correl - mean(muY_for_cross_correl,'omitnan');
                    [r_temp, lags] = xcorr(muX_for_cross_correl,muY_for_cross_correl,...
                        round(fsamp*0.1),"normalized");
                    [max_corr_value, max_corr_lag] = max(r_temp);
                    max_corr_lag = (max_corr_lag - round(length(lags)/2) ) /fsamp ;

                    correl_values_temp_per_trial(triali) = max_corr_value;
                    correl_lags_temp_per_trial(triali) = max_corr_lag;
                end

                flexibility_between_MU_pairs{pair_iter,'cross_corr_trials_mean'} = {mean(correl_values_temp_per_trial,'omitnan')};
                flexibility_between_MU_pairs{pair_iter,'cross_corr_lag_trials_mean'} = {mean(correl_lags_temp_per_trial,'omitnan')};

            end % end of "for each MU Y"
        end % end of "for each MU X"

        close(progressBarMUPair)
        files_each_condition{conditioni}.flexibility_MU_pairs = flexibility_between_MU_pairs;

    end % end of "represented VS unrepresented condition"

    vars_to_save = {"flexibility_between_MU_pairs"...
        };
    save(strcat(output_folder,"/",condition_names{conditioni},"/flexibility_MU_pairs_output.mat"),...
        vars_to_save{:});

    close(progressBarCondition)

end % end of "for each condition"

%% GENERAL TABLE FOR ONE SUBJECT, ALL CONDITIONS, AND SAVE CSV

clearvars flexibility_general_results_table empty_row

colNames_bigtable = colNames; % from the "each mu pair" column names
angle_or_cosine_dist = mu_pair_table.Properties.VariableNames{end};
if isstring(angle_or_cosine_dist)
    angle_or_cosine_dist = convertStringsToChars(angle_or_cosine_dist);
elseif ischar(angle_or_cosine_dist)
end
colNames_bigtable = [{'Subject','Condition'},colNames_bigtable,...
    {'dotprod',angle_or_cosine_dist}];
nb_of_cols = numel(colNames_bigtable);
empty_row = {};
for coli=1:nb_of_cols
    empty_row{coli} = {nan};
end

% Fill a matrix, then convert it to a table (much faster than
% building the table row by row)
flexibility_general_results_temp_mat = [];
% Get the flexibility matrix for each condition first
% subject_idx = already existing because in loaded files

muscle_names_strings = {}; % store muscle names

for conditioni=1:nb_of_conditions
    % first get it as a mat, because much faster than using table
    flexibility_mat_temp = files_each_condition{conditioni}.flexibility_MU_pairs;
    % Every variable is a double, except muscle names which are string.
    % Convert muscle names into NaN values, and re-assign them later
    muscle_names_strings{conditioni} = [flexibility_mat_temp.MU_x_muscle,...
        flexibility_mat_temp.MU_y_muscle,...
        flexibility_mat_temp.MU_pair_muscles];
    flexibility_mat_temp.MU_x_muscle(:) = {nan};
    flexibility_mat_temp.MU_y_muscle(:) = {nan};
    flexibility_mat_temp.MU_pair_muscles(:) = {nan};
    %
    flexibility_mat_temp = table2cell(flexibility_mat_temp);
    flexibility_mat_temp = cell2mat(flexibility_mat_temp);
    nb_rows_temp = size(flexibility_mat_temp,1);
    flexibility_general_results_temp_mat = [flexibility_general_results_temp_mat;...
        [ones(nb_rows_temp,1).*subject_idx,...
        ones(nb_rows_temp,1).*conditioni,...
        flexibility_mat_temp]];
end
% add nan columns for dot prod and angle
flexibility_general_results_temp_mat = [flexibility_general_results_temp_mat,...
    nan(size(flexibility_general_results_temp_mat,1),2)];

% Add the results from the factor analysis
row_nb_for_mu_x_match_idx = 4;
row_nb_for_mu_y_match_idx = 7;
mu_pair_table_as_mat = table2cell(mu_pair_table);
mu_pair_table_as_mat = cell2mat(mu_pair_table_as_mat);
mu_pair_table_as_mat(:,[1,3]) = []; % keep only matched idx (cols 1 and 2) and dot prod and angle (col 3 and 4)
for rowi=1:size(flexibility_general_results_temp_mat,1)
    pair_temp = [...
        flexibility_general_results_temp_mat(rowi,row_nb_for_mu_x_match_idx),...
        flexibility_general_results_temp_mat(rowi,row_nb_for_mu_y_match_idx)];
    if sum(pair_temp==0)<1 % if no unmatched MU (otherwise not possible to know which MU pair it is)
        idx_unit_x = [mu_pair_table_as_mat(:,1)]==pair_temp(1);
        idx_unit_x = idx_unit_x + ...
            ([mu_pair_table_as_mat(:,2)]==pair_temp(1));
        idx_unit_y = [mu_pair_table_as_mat(:,1)]==pair_temp(2);
        idx_unit_y = idx_unit_y + ...
            ([mu_pair_table_as_mat(:,2)]==pair_temp(2));
        idx_unit_x = find(idx_unit_x>=1); % correspond to rows in mu_pair_table
        idx_unit_y = find(idx_unit_y>=1); % correspond to rows in mu_pair_table
        if ~isempty(idx_unit_x) && ~isempty(idx_unit_y)
            temp_mu_pair_row = intersect(idx_unit_x,idx_unit_y);
            if numel(temp_mu_pair_row) > 1
                disp("Conflict");
            end
            if ~isempty(temp_mu_pair_row)
                flexibility_general_results_temp_mat(rowi,nb_of_cols-1) = ... %second-to-last col is for dot prod
                    mu_pair_table_as_mat(temp_mu_pair_row,3);
                flexibility_general_results_temp_mat(rowi,nb_of_cols) = ... %last col is for angle
                    mu_pair_table_as_mat(temp_mu_pair_row,4);
                %                 % if dot prod = 0 and angle = 0, it means that at least one
                %                 % MU from the pai was not matched in the reference
                %                 % condition used for dot prod and angle calculation (it was
                %                 % matched, but not on the condition used for FA)
                %                 if flexibility_general_results_temp_mat(rowi,12) < 0.001 &&...
                %                         flexibility_general_results_temp_mat(rowi,12) > -0.001
                %                     flexibility_general_results_temp_mat(rowi,12) = nan; %12th col is for dot prod
                %                 end
                %                 if flexibility_general_results_temp_mat(rowi,13) < 0.001
                %                     flexibility_general_results_temp_mat(rowi,13) = nan; %13th col is for angle
                %                 end
                continue;
            end
        end
        %     else % else of "if MUs from the pair are matched"
        %         flexibility_general_results_temp_mat(rowi,12) = nan; %12th col is for dot prod
        %         flexibility_general_results_temp_mat(rowi,13) = nan; %13th col is for angle
    end % end of "if MUs from the pair are matched"
end % end of "for each row"

% Save as big table
flexibility_general_results_table = num2cell(flexibility_general_results_temp_mat);
flexibility_general_results_table = cell2table(flexibility_general_results_table);
flexibility_general_results_table.Properties.VariableNames = colNames_bigtable;
% Give back the muscle name strings to the table
full_muscle_names_list = {};
for conditioni=1:nb_of_conditions
    full_muscle_names_list = [full_muscle_names_list;...
        muscle_names_strings{conditioni}];
end
flexibility_general_results_table.MU_x_muscle = full_muscle_names_list(:,1);
flexibility_general_results_table.MU_y_muscle = full_muscle_names_list(:,2);
flexibility_general_results_table.MU_pair_muscles = full_muscle_names_list(:,3);

vars_to_save = {"flexibility_general_results_table";...
    "flexibility_general_results_temp_mat"...
    };
save(strcat(output_folder,"/flexibility_general_table_output.mat"),...
    vars_to_save{:});
writetable(flexibility_general_results_table,strcat(output_folder,...
    "/flexibility_general_table_output.csv"));
