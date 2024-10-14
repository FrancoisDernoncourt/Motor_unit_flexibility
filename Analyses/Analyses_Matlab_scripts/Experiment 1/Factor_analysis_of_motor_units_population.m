clear all
close all
clc
set(0,'DefaultFigureWindowStyle','docked')

%%

subject_idx = 6;
muscle_name = "VL"; %"VL"
title_suffix = "";
condition_names = {"plateau"}; %{"CONT";"PAIN";"WASH"}; %{"CONT";"PAIN"}
corresponding_columns_in_MUs_matched = [4]; %,5,6]; %[4,5]
nb_of_conditions = numel(corresponding_columns_in_MUs_matched);
force_reversed = false;
grids_to_load = 1:6; %1:4; %1:6
grids_to_discard = [1,2]; %[1,2]; %[5,6]; %[];

need_matched_MUs = false;

% Saving options
output_folder = strcat("_Rsquared_per_nb_of_factor_S",num2str(subject_idx));

% FA OPTIONS
max_factors_to_try = 10; %10
FA_on_only_matched_MUs = false;
remove_discontinuous_MUs = true;

% Creating a random dataset
create_shuffled_surrogate = true; % if true, will perform the bootstrap shuffling => both for the identification of nb of factors, and for the clustering
use_downsampling_for_surrogate = true; % recommended, because otherwise very long
shuffling_iterations = 20; %100;

percentile_random_threshold = 95;
mse_threshold_for_straight_curve = 5*(1e-3); % from Cheung 2005 (they used 5*(1e-5))

% Filter discharge rate
% From DelVecchio code neural modules
fsamp = 2048; %set your fsamp;
Wind_s = 0.4;  % hanning window duration. 0.4 for 2.5hz low-pass, 0.2 for 5hz low-pass
HanningW = 2/round(fsamp*Wind_s)*hann(round(fsamp*Wind_s)); %unitary area
detrend_DR_for_each_window = true; % detrend the DR separately for each window
% Basically mostly removes the spike-frequency adaptation
downsample_signal = true;
downsampling_factor = 20;
extend_surrogate_data = false; % if true, will extend surrogate signals to prevent artifacts/drops in DR from the filtering process

resize_concatenated_signals_to_remove_discontinuous_sections = false;
% if true, will shorten the "concatenated" signal to remove all NAN
% values. All MUs thus become continuous
% If false, will use the same approach as usual (entire windows,
% discontinuous MUs not integrated to FA but correlated with factors)
% (DRs look better with second option)

% % % Internal stuff to change colors, titles and folders
if FA_on_only_matched_MUs && need_matched_MUs
    output_folder = strcat(output_folder,"_matched_MUs");
    title_suffix = strcat(title_suffix,"(Matched MUs for FA");
else
    output_folder = strcat(output_folder,"_matched_and_unmatched_MUs");
    title_suffix = strcat(title_suffix,"(Matched and unmatched MUs for FA");
end
if resize_concatenated_signals_to_remove_discontinuous_sections
    output_folder = strcat(output_folder,"_resized_windows_remove_discontinuities");
    title_suffix = strcat(title_suffix," ; windows resized to remove samples with discontinuous DRs");
else
    output_folder = strcat(output_folder,"");
end
title_suffix = strcat(title_suffix,")");
% output_folder = strcat(output_folder,"_clustering_factor_likelihood_to_share_",num2str(arbitrary_clustering_threshold_modulation));

if remove_discontinuous_MUs
    threshold_for_discontinuous_discharge_rate = round(fsamp*Wind_s); % 0.5 second
else
    threshold_for_discontinuous_discharge_rate = inf;
end

%% GET PATHS AND FILENAMES

paths_string_to_load = {};
files_string_to_load = {};

% signal files
for filei=1:nb_of_conditions
    [files_string_to_load{end+1}, paths_string_to_load{end+1}] = ...
        uigetfile(".mat",strcat("Please load the signal file corresponding to condition #",num2str(filei)), ...
        strcat("Please load the signal file corresponding to condition #",num2str(filei)), ...
        "Multiselect","off");
end

% window limitis files
for filei=1:nb_of_conditions
    [files_string_to_load{end+1}, paths_string_to_load{end+1}] = ...
        uigetfile(".mat",strcat("Please load the file with window edges corresponding to condition #",num2str(filei)), ...
        strcat("Please load the signal file window edges corresponding to condition #",num2str(filei)), ...
        "Multiselect","off");
end

% matched file
if need_matched_MUs
    [files_string_to_load{end+1}, paths_string_to_load{end+1}] = ...
        uigetfile(".mat",strcat("Please load the matched MUs file"), ...
        strcat("Please load the matched MUs file"), ...
        "Multiselect","off");
end

%% ACTUALLY LOAD FILES
mkdir(output_folder)
load(strcat(paths_string_to_load{end},files_string_to_load{end})); % Load MUs matched file

% Remove grids to discard from the beginning
if ~isempty(grids_to_discard)
    for gridi_discard = 1:numel(grids_to_discard)
        grids_to_load(find(grids_to_load==grids_to_discard(gridi_discard))) = [];
    end
end

for filei=1:numel(files_string_to_load)
    if filei <= nb_of_conditions % signal files
        load(strcat(paths_string_to_load{filei},files_string_to_load{filei}));
        % Remove empty grids from grids to load
        grids_to_load_for_file = grids_to_load;
        for gridi=1:numel(grids_to_load)
            if isempty(edition.Pulsetrain{grids_to_load(gridi)})
                grids_to_load_for_file(find(grids_to_load_for_file==grids_to_load(gridi))) = [];
            end
        end

        %files_each_condition{filei}.force = signal.path;
        files_each_condition{filei}.force = signal.data(signal.ngrid*64+1,:);
        if force_reversed
            files_each_condition{filei}.force = files_each_condition{filei}.force .* (-1);
        end
        files_each_condition{filei}.time = (1:length(files_each_condition{filei}.force))./fsamp;
        files_each_condition{filei}.binary_spike_trains = edition.Dischargetimes(grids_to_load_for_file,:);
        if need_matched_MUs
            files_each_condition{filei}.condition = MUs_matched.Properties.VariableNames{corresponding_columns_in_MUs_matched(filei)};
        end

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

%% IF DOWNSAMPLING, CHANGE DIRECTLY IN THE "FILES TO LOAD" AND ADAPT EVERY OTHER VALUE ACCORDINGLY
if downsample_signal
    for filei=1:numel(files_each_condition)
        files_each_condition{filei}.force =  downsample(files_each_condition{filei}.force, downsampling_factor);
        files_each_condition{filei}.time =  downsample(files_each_condition{filei}.time, downsampling_factor);
        for mui=1:numel(files_each_condition{filei}.binary_spike_trains)
            if ~isempty(files_each_condition{filei}.binary_spike_trains{mui})
                files_each_condition{filei}.binary_spike_trains{mui} = round( files_each_condition{filei}.binary_spike_trains{mui} ./ downsampling_factor);
            end
        end
        window_lims_files_each_condition{filei}.limofEachWindow = round( window_lims_files_each_condition{filei}.limofEachWindow ./ downsampling_factor);
        window_lims_files_each_condition{filei}.limofEachWindow_continuous(:,1) = ceil(window_lims_files_each_condition{filei}.limofEachWindow_continuous(:,1)./downsampling_factor);
        window_lims_files_each_condition{filei}.limofEachWindow_continuous(:,2) = floor(window_lims_files_each_condition{filei}.limofEachWindow_continuous(:,2)./downsampling_factor);
    end
    limOfEachWindow = round(limOfEachWindow./downsampling_factor);
    limofEachWindow_continuous(:,1) = ceil(limofEachWindow_continuous(:,1)./downsampling_factor);
    limofEachWindow_continuous(:,2) = floor(limofEachWindow_continuous(:,2)./downsampling_factor);

    %
    fsamp = fsamp / downsampling_factor;
    threshold_for_discontinuous_discharge_rate = round(threshold_for_discontinuous_discharge_rate / downsampling_factor);
    HanningW = 2/round(fsamp*Wind_s)*hann(round(fsamp*Wind_s)); %unitary area
    HanningW = HanningW ./sum(HanningW); % making sure the sum = 1
end

%% CREATE ACTUAL SIGNAL FOR EACH FILE

if need_matched_MUs
    % change empty cells to zeros in MUs_matched
    MUs_matched = MUs_matched(1:max([MUs_matched{:,"MU_idx"}{:}]),:);
    for vari = 1:size(MUs_matched,2)
        MUs_matched{find(cellfun(@isempty, MUs_matched{:,vari})),vari} = {0};
    end
end

nb_of_grids = numel(grids_to_load_for_file);
for conditioni=1:nb_of_conditions
    binary_matrix = [];
    smoothed_DR = [];
    smoothed_DR_discontinuous_nans = [];
    mu_total = 0;
    corresponding_MUs_for_condition = [];
    for gridi=1:nb_of_grids
        grid_idx = grids_to_load_for_file(gridi);
%         if sum(grid_idx==grids_to_discard) % skip grids to ignore = but
%         should be already removed anyway from the loading part of the
%         script
%             continue
%         end
        nonempty_mu_idx_in_grid = find( ~cellfun(@isempty, files_each_condition{conditioni}.binary_spike_trains(gridi,:) ) );
        % At least 10 spikes are necessary (prevent problems for FA)
        mu_too_few_spikes = [];
        for mui=1:numel(nonempty_mu_idx_in_grid)
            if numel(files_each_condition{conditioni}.binary_spike_trains{gridi,mui}) <= 10
                mu_too_few_spikes(end+1) = mui;
            end
        end
        nonempty_mu_idx_in_grid(mu_too_few_spikes) = [];
        clearvars mu_too_few_spikes
        %
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
            if need_matched_MUs
                corresponding_MUs_matched_row = ...
                    find( [MUs_matched{:,corresponding_columns_in_MUs_matched(conditioni)}{:}] == temp_idx);
                if ~isempty(corresponding_MUs_matched_row)
                    corresponding_MUs_for_condition(mu_total,2) = corresponding_MUs_matched_row;
                else
                    corresponding_MUs_for_condition(mu_total,2) = 0;
                end
            end
        end % end of "for each MU"
    end % end of "for eadh grid"

    smoothed_DR_concatenated_with_nans = [];
    smoothed_DR_concatenated_no_nans = [];
    binary_matrix_concatenated = [];
    limofEachWindow = window_lims_files_each_condition{conditioni}.limofEachWindow;
    limofEachWindow_continuous = window_lims_files_each_condition{conditioni}.limofEachWindow_continuous;
    nb_of_windows = size(limofEachWindow,1);
    for windowi=1:nb_of_windows
        current_length = size(binary_matrix_concatenated,2);
        new_length = length(limofEachWindow(windowi,1):limofEachWindow(windowi,2));
        binary_matrix_concatenated(:, ...
            current_length+1:current_length+new_length) = ...
            binary_matrix(:,limofEachWindow(windowi,1):limofEachWindow(windowi,2));
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
     clearvars smoothed_DR_current_window

     if resize_concatenated_signals_to_remove_discontinuous_sections % If true, modify window limits (window resize)
         % Remove NaN idx from the concatenated signal (in a copy of the
         % signal)
         % (only in "smoothed_DR_concatenated_with_nans")
         nan_idx = [];
         for mui = 1:size(smoothed_DR_concatenated_with_nans,1)
             nan_idx(mui,:) = isnan(smoothed_DR_concatenated_with_nans(mui,:));
         end
         nan_idx_full = sum(nan_idx);
         samples_to_remove = find(nan_idx_full>=1);
         smoothed_DR_concatenated_with_nans(:,samples_to_remove) = [];


         binary_matrix_concatenated(:,samples_to_remove) = [];
         nb_removed_idx_each_window = [];
         limofEachWindow_continuous_temp = limofEachWindow_continuous;
         for windowi=1:nb_of_windows
             temp_samples_to_remove = samples_to_remove(samples_to_remove >= limofEachWindow_continuous_temp(windowi,1));
             temp_samples_to_remove = temp_samples_to_remove(temp_samples_to_remove < limofEachWindow_continuous_temp(windowi,2));
             nb_removed_idx_each_window(windowi) = numel(temp_samples_to_remove);
             limofEachWindow_continuous(windowi+1:end,1) = limofEachWindow_continuous(windowi+1:end,1)-nb_removed_idx_each_window(windowi);
             limofEachWindow_continuous(windowi:end,2) = limofEachWindow_continuous(windowi:end,2)-nb_removed_idx_each_window(windowi);
         end
         files_each_condition{conditioni}.nb_removed_idx_each_window = nb_removed_idx_each_window;
         clearvars limofEachWindow_continuous_temp temp_samples_to_remove nan_idx samples_to_remove
     end
     %

    corresponding_MUs_for_condition(:,3) = NaN; % Put "status" of MU in the third column
    % (display matched MU idx, not idx based on RT)
    % 1 = matched and continuous MUs = used for FA
    % 0 = unmatched and continuous MUs = not used for FA but correlated with latents
    % -1 = matched and discontinuous MUs = not used for FA but correlated with latents
    % -2 = unmatched and discontinuous MUs = not used for FA but correlated with latents
    % MUs from discarded grids = should not be present in the MU list
    % anyway
    for mui=1:size(smoothed_DR_concatenated_with_nans,1)
        if sum(isnan(smoothed_DR_concatenated_with_nans(mui,:))) < threshold_for_discontinuous_discharge_rate
            if corresponding_MUs_for_condition(mui,2) >= 1 % matched and continuous MUs = used for FA
                corresponding_MUs_for_condition(mui,3) = 1;
            else % unmatched and continuous MUs = not used for FA but correlated with latents
                corresponding_MUs_for_condition(mui,3) = 0;
            end
        else
            if corresponding_MUs_for_condition(mui,2) >= 1 % matched and discontinuous MUs = not used for FA but correlated with latents
                corresponding_MUs_for_condition(mui,3) = -1;
            else % unmatched and discontinuous MUs = not used for FA but correlated with latents
                corresponding_MUs_for_condition(mui,3) = -2;
            end
        end
    end

    files_each_condition{conditioni}.window_edges = limofEachWindow;
    files_each_condition{conditioni}.window_edges_continuous = limofEachWindow_continuous;
    files_each_condition{conditioni}.binary_matrix = binary_matrix;
    files_each_condition{conditioni}.binary_matrix_concatenated = binary_matrix_concatenated;
    files_each_condition{conditioni}.smoothed_DR = smoothed_DR;
    files_each_condition{conditioni}.smoothed_DR_concatenated_no_nans = smoothed_DR_concatenated_no_nans;
    files_each_condition{conditioni}.smoothed_DR_discontinuous_nans = smoothed_DR_discontinuous_nans;
    files_each_condition{conditioni}.smoothed_DR_concatenated_with_nans = smoothed_DR_concatenated_with_nans;

    % Convert "corresponding MUs" into table to make it easier to work with it
    columnNames = {'MU_grid_and_idx_in_file', 'MU_idx_in_matched_MUs', 'Status'};
    corresponding_MUs_for_condition = array2table(corresponding_MUs_for_condition, 'VariableNames', columnNames);
    col_status_temp = num2cell(corresponding_MUs_for_condition.Status);
    for mui = 1:size(col_status_temp,1)
        switch col_status_temp{mui}
            case 1
                col_status_temp{mui} = "matched_continuous";
            case 0
                col_status_temp{mui} = "unmatched_continuous";
            case -1
                col_status_temp{mui} = "matched_discontinuous";
            case -2
                col_status_temp{mui} = "unmatched_discontinuous";
            otherwise
                continue
        end
    end
    corresponding_MUs_for_condition.Status = col_status_temp;
    files_each_condition{conditioni}.corresponding_MUs_for_condition = corresponding_MUs_for_condition;

%     % Sanity check - commented everything to redisplay it later, to give MU
%     % color according to correlation with factors
%     close(figure(conditioni))
%     figure(conditioni)
%     color_for_each_MU_temp = turbo(size(smoothed_DR_concatenated_with_nans,1));
%     color_for_each_MU_temp(:,4) = 1;
%     for mui=1:size(smoothed_DR_concatenated_with_nans,1)
% %         if resize_concatenated_signals_to_remove_discontinuous_sections
% %             plot( (1:size(smoothed_DR_concatenated_with_nans,2))./fsamp,...
% %                 smoothed_DR_concatenated_with_nans(mui,:) );
% %         else
% %             plot( (1:size(smoothed_DR_concatenated_no_nans,2))./fsamp,...
% %                 smoothed_DR_concatenated_no_nans(mui,:) );
% %         end
%         if sum(isnan(smoothed_DR_concatenated_with_nans(mui,:))) >= 1
%             color_for_each_MU_temp(mui,:) = [.25,.25,.25,0.3];
%         end
%         plot( (1:size(smoothed_DR_concatenated_with_nans,2))./fsamp,...
%             smoothed_DR_concatenated_with_nans(mui,:), 'Color', color_for_each_MU_temp(mui,:) );
%         hold on
%     end
%     for windowi=1:nb_of_windows
%         line([limofEachWindow_continuous(windowi,2),limofEachWindow_continuous(windowi,2)]./fsamp, ...
%             ylim,'LineStyle',':','Color','Black','LineWidth',3);
%     end
%     ylabel("Discharge rate (spike/s)");
%     xlabel("Time (s)");
%     title({strcat(subject_idx," - condition ", num2str(conditioni), " (", condition_names{conditioni},") - DR sanity check (only continuous MUs)"); ...
%         title_suffix;...
%         strcat("Discontinuous MUs (not used for FA) in grey ; threshold for discontinuity = ",...
%         num2str(round(threshold_for_discontinuous_discharge_rate*10/fsamp)/10),"s")});
%     box off
%     saveas(figure(conditioni*1),strcat(output_folder,"/DRs_",condition_names{conditioni},".png"));

end % enf of "for each condition"

%% CREATE SURROGATE (RANDOMLY SHUFFLED) DATA FOR CLUSTERING AND FACTOR SELECTION
threshold_for_ISI_too_short = 0.02; % Check for ISIs that are very small (can happen because of concatenation)
% threshold in seconds. 0.02 is 50hz IDR
threshold_for_ISI_too_short = threshold_for_ISI_too_short * fsamp;

%only for 1st file, because won't be used for more than 1 file anyway
surrogate_data = struct();
surrogate_data.binary_matrix = cell(shuffling_iterations,1);
surrogate_data.smoothed_DRs = cell(shuffling_iterations,1);

waitbar_text = "Generating surrogate/baseline dataset (randomly shuffled spike trains)";
wait_bar = waitbar(0,waitbar_text);

for shuffli=1:shuffling_iterations
    for mui=1:size(files_each_condition{1}.binary_matrix_concatenated,1)
        % binary mat
        surrogate_data.binary_matrix{shuffli}(mui,:) = files_each_condition{1}.binary_matrix_concatenated(mui,:);
        cumsum_binary = cumsum(surrogate_data.binary_matrix{shuffli}(mui,:));
        unique_idx_temp = unique(cumsum_binary);
        ISIs = [];
        for spikei=1:length(unique_idx_temp)-1
            ISIs(spikei) = sum(cumsum_binary==unique_idx_temp(spikei));
        end
        % Check for ISIs that are very small (can happen because of
        % concatenation), and remove those spikes
        spikes_to_remove_ISI_too_short = find(ISIs<threshold_for_ISI_too_short);
        if ~isempty(spikes_to_remove_ISI_too_short)
            ISIs(spikes_to_remove_ISI_too_short) = [];
        end
        ISIs_shuffled = ISIs(randperm(length(ISIs))); % Random shuffling
        cumsum_ISIs_shuffled = cumsum(ISIs_shuffled);
        surrogate_data_temp = zeros(1,length(surrogate_data.binary_matrix{shuffli}(mui,:)));
        for spikei=1:length(cumsum_ISIs_shuffled)
            surrogate_data_temp(cumsum_ISIs_shuffled(spikei)) = 1;
        end
        surrogate_data.binary_matrix{shuffli}(mui,:) = surrogate_data_temp;

        % Extend with random edges to prevent artifacts from filtering
        if extend_surrogate_data
            edge_extension = 2; % in seconds
            random_edge = {};
            % 1st edge (beginning)
            random_edge_start = randi(length(surrogate_data_temp) - ceil(fsamp*edge_extension));
            random_edge{1} = surrogate_data_temp(random_edge_start+1:random_edge_start+ceil(fsamp*edge_extension));
            %             random_edge{1}([ (1:round(threshold_for_ISI_too_short)) , ...
            %                 (end-round(threshold_for_ISI_too_short):end) ] ) = 0;
            random_edge{1}([ (1:median(ISIs)) , ...
                ( length(random_edge{1})-floor(median(ISIs)) : length(random_edge{1}) ) ] ) = 0;
            % 2nd edge (end)
            random_edge_start = randi(length(surrogate_data_temp) - ceil(fsamp*edge_extension));
            random_edge{2} = surrogate_data_temp(random_edge_start+1:random_edge_start+ceil(fsamp*edge_extension));
            %             random_edge{2}([ (1:round(threshold_for_ISI_too_short)) , ...
            %                 (end-round(threshold_for_ISI_too_short):end) ] ) = 0;
            random_edge{2}([ (1:median(ISIs)) , ...
                ( length(random_edge{2})-floor(median(ISIs)) : length(random_edge{2}) ) ] ) = 0;
            % Add edges
            surrogate_data_temp(1) = 1; % Make the fir
            surrogate_data_temp(end+1:end+length(random_edge{2})) = random_edge{2};
            surrogate_data_temp(1:end+length(random_edge{1})) = [random_edge{1},surrogate_data_temp];
            % Add a first "true" firing around the beginning of the surrogate here,
            % because the shuffled ISI assumes that there would be a spike
            % there (otherwise all surrogate MUs have a drop in DR here)
            add_idx = round(length(random_edge{1}) + (randn(1)*min(ISIs)));
            surrogate_data_temp(add_idx) = randi([0,1])*0.5; %randi([0,1]);
            add_idx = round(length(random_edge{2}) + (randn(1)*min(ISIs)));
            surrogate_data_temp(end-add_idx) = randi([0,1])*0.5; %randi([0,1]);
        else
            random_edge{1} = 1:floor(fsamp/2);
            random_edge{2} = 1:floor(fsamp/2);
        end

        % Filter = smoothed DR
        smoothed_surrogate_temp = filtfilt(HanningW,1,surrogate_data_temp*fsamp);
        % remove artifacted edges
        smoothed_surrogate_temp([ (1:length(random_edge{1})) , ...
            (end-length(random_edge{2}):end) ]) = [];
        % Duplicate or remove a sample to have the same size
        if length(smoothed_surrogate_temp) < length(surrogate_data.binary_matrix{shuffli}(mui,:))
            duplicate_sample = randi([1,length(smoothed_surrogate_temp)]);
            smoothed_surrogate_temp(duplicate_sample+1:end+1) = smoothed_surrogate_temp(duplicate_sample:end);
        elseif length(smoothed_surrogate_temp) > length(surrogate_data.binary_matrix{shuffli}(mui,:))
            smoothed_surrogate_temp(1) = [];
        end

        % Add NANs to discontinuous gaps
        for spiki=1:length(ISIs_shuffled)
            if ISIs_shuffled(spiki) > threshold_for_discontinuous_discharge_rate
                if spiki == 1 % if first spike
                    smoothed_surrogate_temp(1:cumsum_ISIs_shuffled(spiki)) = nan;
                elseif spiki == length(cumsum_ISIs_shuffled) % if last spike
                    smoothed_surrogate_temp(cumsum_ISIs_shuffled(spiki-1):end) = nan;
                else
                    smoothed_surrogate_temp(cumsum_ISIs_shuffled(spiki-1):cumsum_ISIs_shuffled(spiki)) = nan;
                end
            end
        end
        %

        % Detrend
        if detrend_DR_for_each_window
            smoothed_surrogate_temp = ...
                detrend(smoothed_surrogate_temp,'omitnan');
        end
        %

        % Assign the variable
        if mui == 1 % just making sure the signal pasted has the right length
            surrogate_data.smoothed_DRs{shuffli}(mui,:) = smoothed_surrogate_temp; 
        else
            if length(smoothed_surrogate_temp) < length(surrogate_data.smoothed_DRs{shuffli})
                % pad with zeros if missing a few samples
                smoothed_surrogate_temp( ...
                    end+1 : length(surrogate_data.smoothed_DRs{shuffli}) ) = 0;
            end
            surrogate_data.smoothed_DRs{shuffli}(mui,:) = 0;
            surrogate_data.smoothed_DRs{shuffli}(mui,:) = smoothed_surrogate_temp(...
                1:length(surrogate_data.smoothed_DRs{shuffli}) );
        end
    end % end of "for each MU"

    % plot(surrogate_data.smoothed_DRs{shuffli}');
    waitbar(shuffli/shuffling_iterations, wait_bar, ...
        waitbar_text);
end % end of "for each shuffle"

close(wait_bar)

%% END OF PREPROCESSING

for conditioni=1:nb_of_conditions

    % Remove discontinuous and unmatched MUs
    if resize_concatenated_signals_to_remove_discontinuous_sections
        temp_smoothed_DR_for_FA = files_each_condition{conditioni}.smoothed_DR_concatenated_with_nans; % use shortened version with NANs removed
    else
        temp_smoothed_DR_for_FA = files_each_condition{conditioni}.smoothed_DR_concatenated_no_nans; % use full version (but still no NANs because no replacement of discontinuous times with NANs)
    end

    idx_remove_for_FA = [];
    for mui=size(files_each_condition{conditioni}.corresponding_MUs_for_condition,1):-1:1 % going in reverse because removing
        if files_each_condition{conditioni}.corresponding_MUs_for_condition.Status{mui} ~= "matched_continuous"
            if ~FA_on_only_matched_MUs &&... % if FA on unmatched MUs also
                    files_each_condition{conditioni}.corresponding_MUs_for_condition.Status{mui} ~= "unmatched_continuous"
                temp_smoothed_DR_for_FA(mui,:) = [];
                idx_remove_for_FA(end+1) = mui;
                continue
            end
        end
        % temp_smoothed_DR_for_FA(mui,isnan(temp_smoothed_DR_for_FA(mui,:))) = 0; % remove NAN values
    end
    temp_smoothed_DR_for_FA = temp_smoothed_DR_for_FA';
   

    %% GET R² OF RECONSTRUCTED VS ORIGINAL DATA - USING ACTUAL FA
    % ALSO GET THE DROP IN INDEPENDENT VARIANCE
    surrogate_rsquared = zeros(shuffling_iterations,max_factors_to_try+1);
    surrogate_independent_var = ones(shuffling_iterations,max_factors_to_try+1);
    original_rsquared = zeros(1,max_factors_to_try+1);
    original_independent_var = ones(1,max_factors_to_try+1); 

    waitbar_text = "Factor analysis for original and surrogate data, for different numbers of factors";
    wait_bar = waitbar(0, waitbar_text);

    total_iter_FA(1) = 0;
    total_iter_FA(2) = max_factors_to_try * (shuffling_iterations + 1);

    % Original data
    for latent_factori=1:max_factors_to_try
        total_iter_FA(1) = total_iter_FA(1) + 1;

        [factor_weights, specific_vars , ~, ~, factor_time_vec] = ...
            factoran(temp_smoothed_DR_for_FA,latent_factori,"xtype","data","rotate","promax","maxit",1e5);
        reconstructed_data_temp = factor_weights * factor_time_vec';
        reconstructed_data_temp = reconstructed_data_temp';
        per_mu_rsquared_temp = [];
        for mui=1:size(temp_smoothed_DR_for_FA,2)
            per_mu_rsquared_temp(mui) = corr( ...
                temp_smoothed_DR_for_FA(:,mui), ...
                reconstructed_data_temp(:,mui)).^2;
        end
        original_rsquared(1,latent_factori+1) = mean(per_mu_rsquared_temp);
        original_independent_var(1,latent_factori+1) = mean(specific_vars);
        waitbar(total_iter_FA(1)/total_iter_FA(2), wait_bar, ...
            waitbar_text);
    end

    % Surrogate data (a bit long... because one FA per surrogate dataset
    % (one per shuffle) x one FA per number of latent factors to include)
    for shuffli=1:shuffling_iterations
        surrogate_for_FA_temp = surrogate_data.smoothed_DRs{shuffli};
%         surrogate_nans_to_delete = isnan(surrogate_for_FA_temp);
%         surrogate_nans_to_delete = find(sum(surrogate_nans_to_delete)>1);
        surrogate_for_FA_temp(isnan(surrogate_for_FA_temp)) = 0; % Just a quick fix for some situations where some NANs occur
        surrogate_for_FA_temp(idx_remove_for_FA,:) = [];
        surrogate_for_FA_temp = surrogate_for_FA_temp';
        for latent_factori=1:max_factors_to_try
            total_iter_FA(1) = total_iter_FA(1) + 1;

            [factor_weights, specific_vars , ~, ~, factor_time_vec] = ...
                factoran(surrogate_for_FA_temp,latent_factori,"xtype","data","rotate","promax","maxit",1e5);
            reconstructed_data_temp = factor_weights * factor_time_vec';
            reconstructed_data_temp = reconstructed_data_temp';
            per_mu_rsquared_temp = [];
            for mui=1:size(surrogate_for_FA_temp,2)
                per_mu_rsquared_temp(mui) = corr( ...
                    surrogate_for_FA_temp(:,mui), ...
                    reconstructed_data_temp(:,mui)).^2;
            end
            surrogate_rsquared(shuffli,latent_factori+1) = mean(per_mu_rsquared_temp);
            surrogate_independent_var(shuffli,latent_factori+1) = mean(specific_vars);
            waitbar(total_iter_FA(1)/total_iter_FA(2), wait_bar, ...
                waitbar_text);
        end
    end
    %
    surrogate_rsquared_mean = mean(surrogate_rsquared);
    surrogate_independent_var_mean = mean(surrogate_independent_var);

    close(wait_bar)

    % % % % %
    %% Different mehtods to find the appropriate number of factors to use
    % % % % %

    % linear fit error
    linear_fit_mse_of_Rsquared_curve = [];
    slope_of_Rsquared_curve = [];
    for i = 2:length(original_rsquared)-1
        x_temp = [i-1:i+1];
        y_temp = original_rsquared(i-1:i+1);
        fit_coefs_temp = polyfit(x_temp, y_temp, 1);
        slope_of_Rsquared_curve(i-1) = fit_coefs_temp(1);
        fit_prediction_temp = polyval(fit_coefs_temp, x_temp);
        mse_of_fit_temp = mean(sqrt((y_temp - fit_prediction_temp) .^ 2));
        linear_fit_mse_of_Rsquared_curve(i-1) = mse_of_fit_temp;
    end
    linear_fit_mse_of_surrogate_Rsquared_curve = [];
    %slope_of_surrogate_Rsquared = [];
    for i = 2:length(surrogate_rsquared_mean)-1
        x_temp = [i-1:i+1];
        y_temp = surrogate_rsquared_mean(i-1:i+1);
        fit_coefs_temp = polyfit(x_temp, y_temp, 1);
        %slope_of_surrogate_Rsquared(i-1) = fit_coefs_temp(1);
        fit_prediction_temp = polyval(fit_coefs_temp, x_temp);
        mse_of_fit_temp = mean(sqrt((y_temp - fit_prediction_temp) .^ 2));
        linear_fit_mse_of_surrogate_Rsquared_curve(i-1) = mse_of_fit_temp;
    end
    slope_temp = polyfit(1:max_factors_to_try, surrogate_rsquared_mean(1:max_factors_to_try), 1);
    slope_of_surrogate_Rsquared = slope_temp(1);

    % Thresholds options for selecting nb of latents
    %

%     % ::: R² for additional latent for experimental data < 95th percentile % R² for one latent for surrogate data (me)
%     % 1st value is threshold, second value is nb of latents selected by
%     % this threshold
%     threshold_for_nb_of_latent_selection.additional_Rsquared_inferior_to_expected_by_chance{1} = ...
%             prctile(surrogate_rsquared_mean(:,2),percentile_random_threshold); %2 because starts at 2
%     temp_nb_of_latent_selection = ...
%         diff(original_rsquared) < threshold_for_nb_of_latent_selection.additional_Rsquared_inferior_to_expected_by_chance{1};
%     temp_nb_of_latent_selection = min(find(temp_nb_of_latent_selection))-1; % last index to meet the criteria (threshold)
%     if isempty(temp_nb_of_latent_selection)
%         temp_nb_of_latent_selection = 1;
%     end
%     threshold_for_nb_of_latent_selection.additional_Rsquared_inferior_to_expected_by_chance{2} = temp_nb_of_latent_selection;

    % ::: R² > 80% or 90% (Del Vecchio 2023, Torres-Oviedo 2006)
    % = doesn't work with VL
    
    % ::: Additional R² < 5% (Clark 2010)
    threshold_for_nb_of_latent_selection.additional_Rsquared_inferior_to_arbitrary_thresh{1} = ...
        5*0.01;
    temp_nb_of_latent_selection = ...
        diff(original_rsquared) > threshold_for_nb_of_latent_selection.additional_Rsquared_inferior_to_arbitrary_thresh{1};
    temp_nb_of_latent_selection = max(find(temp_nb_of_latent_selection)); % last index to meet the criteria (threshold)
    threshold_for_nb_of_latent_selection.additional_Rsquared_inferior_to_arbitrary_thresh{2} = temp_nb_of_latent_selection;

    % ::: MSE of linear fit < 5 x 10^(-5) (Cheung 2005 = 1.e-5 originally)
    threshold_for_nb_of_latent_selection.mse_of_linear_fit_mse_Rsquared{1} = ...
        mse_threshold_for_straight_curve;
    temp_nb_of_latent_selection = ...
        linear_fit_mse_of_Rsquared_curve < threshold_for_nb_of_latent_selection.mse_of_linear_fit_mse_Rsquared{1};
    temp_nb_of_latent_selection = min(find(temp_nb_of_latent_selection)) - 1; % last index to be above the threshold (hence the -1)
    threshold_for_nb_of_latent_selection.mse_of_linear_fit_mse_Rsquared{2} = temp_nb_of_latent_selection;
    if isempty(temp_nb_of_latent_selection)
        threshold_for_nb_of_latent_selection.mse_of_linear_fit_mse_Rsquared{2} = 99;
    end

    % ::: original-slope drops below 75% of the surrogate-slope for
    % reconstructed data (using R² of reconstructed VS original data)
    % (Cheung 2009 PNAS)
    threshold_for_nb_of_latent_selection.slope_is_below_surrogate_Rsquared{1} = ...
        slope_of_surrogate_Rsquared;
    % Possible to change the threshold...
    temp_nb_of_latent_selection = ...
        slope_of_Rsquared_curve < slope_of_surrogate_Rsquared;
    temp_nb_of_latent_selection = min(find(temp_nb_of_latent_selection)) - 1; % last index to be above the threshold (hence the -1)
    if isempty(temp_nb_of_latent_selection)
        temp_nb_of_latent_selection = 1;
    end
    threshold_for_nb_of_latent_selection.slope_is_below_surrogate_Rsquared{2} = temp_nb_of_latent_selection;

    %%  Plot figure with results
    close(figure(2))
    figure(2)

    default_colors = lines(10);
    color_temp_surrogate = default_colors(2,:);
    temp_marker_color_surrogate = default_colors(2,:)./2;
    color_temp = default_colors(1,:);
    temp_marker_color = default_colors(1,:)./2;

    upper_bound = [prctile(surrogate_rsquared,95)];
    lower_bound = [prctile(surrogate_rsquared,5)];
    % Create a matrix for the fill function
    x = [0:max_factors_to_try];
    fill_x = [x, fliplr(x)];  % fliplr is used to reverse the order for the lower bound
    fill_y = [upper_bound(1:max_factors_to_try+1), fliplr(lower_bound(1:max_factors_to_try+1))];
    % Display the shaded area
    hold on;
    fill(fill_x, fill_y, color_temp_surrogate, 'FaceAlpha', 0.3, 'LineStyle','none');
    plot(0:max_factors_to_try, ...
        surrogate_rsquared_mean(1:max_factors_to_try+1),'--o','LineWidth',4, 'MarkerSize', 10, ...
        'MarkerFaceColor',temp_marker_color_surrogate,...
        'color',color_temp_surrogate);
    % Original data
    plot(0:max_factors_to_try, ...
        original_rsquared(1:max_factors_to_try+1),'-o','LineWidth',4, 'MarkerSize', 10, ...
        'MarkerFaceColor',temp_marker_color,...
        'color',color_temp);
    box off
    xlabel("Number of components")
    ylabel("R² (Coefficient of determination between reconstructed data and true data)")
    ylim([0,1]);

    % Show number of latents to select according to the different methods
%     line([threshold_for_nb_of_latent_selection.additional_Rsquared_inferior_to_expected_by_chance{2},...
%         threshold_for_nb_of_latent_selection.additional_Rsquared_inferior_to_expected_by_chance{2}], ...
%         ylim(), 'LineStyle', ':', 'Color', [0.5, 0, 0.6],'LineWidth', 3);
%     text(threshold_for_nb_of_latent_selection.additional_Rsquared_inferior_to_expected_by_chance{2} + 0.05,...
%         98.5/100, ...
%         strcat('Additional R² < R² for 1 component with surrogate data (', num2str(percentile_random_threshold),' percentile)'));

    line([threshold_for_nb_of_latent_selection.additional_Rsquared_inferior_to_arbitrary_thresh{2},...
        threshold_for_nb_of_latent_selection.additional_Rsquared_inferior_to_arbitrary_thresh{2}], ...
        ylim(), 'LineStyle', ':', 'Color', [0.5, 0, 0.6],'LineWidth', 3);
    text(threshold_for_nb_of_latent_selection.additional_Rsquared_inferior_to_arbitrary_thresh{2} + 0.05,...
        97/100, 'Additional R² < arbitrary threshold (0.05) (Clark 2010)');

    if threshold_for_nb_of_latent_selection.mse_of_linear_fit_mse_Rsquared{2} > max_factors_to_try
        text(max_factors_to_try,0.9, ...
            {'Error (MSE) of linear fit (R² curve) < small threshold (0.005) (curve becomes straight) (Cheung 2005)', ...
            strcat("Nb of factors with this method is > ",num2str(max_factors_to_try)) }, ...
            'HorizontalAlignment','right');
    else
        line([threshold_for_nb_of_latent_selection.mse_of_linear_fit_mse_Rsquared{2},...
            threshold_for_nb_of_latent_selection.mse_of_linear_fit_mse_Rsquared{2}], ...
            ylim(), 'LineStyle', ':', 'Color', [0.5, 0, 0.6],'LineWidth', 3);
        text(threshold_for_nb_of_latent_selection.mse_of_linear_fit_mse_Rsquared{2} + 0.05,...
            95.5/100, 'Error (MSE) of linear fit (R² curve) < small threshold (0.005) (curve becomes straight) (Cheung 2005)');
    end

    line([threshold_for_nb_of_latent_selection.slope_is_below_surrogate_Rsquared{2},...
        threshold_for_nb_of_latent_selection.slope_is_below_surrogate_Rsquared{2}], ...
        ylim(), 'LineStyle', ':', 'Color', [0.5, 0, 0.6],'LineWidth', 3);
    text(threshold_for_nb_of_latent_selection.slope_is_below_surrogate_Rsquared{2} + 0.05,...
        94/100, 'Slope (R² curve) is below surrogate mean slope (Cheung 2009)');

    % Add title and save
    title(strcat("S",num2str(subject_idx), " - R² (reconstructed VS true data) and number of components to select"))
    saveas(gcf(),strcat(output_folder,"/Rsquared_factor_nb_selection_",condition_names{conditioni},".png"))

    %% ORGANIZE AND SAVE OUTPUT

    output_table = table();
    table_vars = {'Subject';'Muscle';...
        'Nb_factors_linear_fit';'Rsquared_associated_linear_fit';...
        'Nb_factors_below_surrogate';'Rsquared_associated_below_surrogate';...
        'Nb_factors_5percent_threshold';'Rsquared_associated_5percent_threshold'};
    empty_row = cell(1,numel(table_vars));
    output_table = cell2table(empty_row,"VariableNames",table_vars);

    output_table.Subject{1} = subject_idx;
    output_table.Muscle{1} = muscle_name;
    output_table.Nb_factors_linear_fit{1} = threshold_for_nb_of_latent_selection.mse_of_linear_fit_mse_Rsquared{2};
    if threshold_for_nb_of_latent_selection.mse_of_linear_fit_mse_Rsquared{2} > max_factors_to_try
        output_table.Rsquared_associated_linear_fit{1} = nan;
    else
        output_table.Rsquared_associated_linear_fit{1} = original_rsquared(...
            output_table.Nb_factors_linear_fit{1} + 1); %+1 because starts at 0
    end
    output_table.Nb_factors_below_surrogate{1} = threshold_for_nb_of_latent_selection.slope_is_below_surrogate_Rsquared{2};
    output_table.Rsquared_associated_below_surrogate{1} = original_rsquared(...
        output_table.Nb_factors_below_surrogate{1}+1); %+1 because starts at 0
    output_table.Nb_factors_5percent_threshold{1} = ...
        threshold_for_nb_of_latent_selection.additional_Rsquared_inferior_to_arbitrary_thresh{2};
    output_table.Rsquared_associated_5percent_threshold{1} = original_rsquared(...
        output_table.Nb_factors_5percent_threshold{1}+1); %+1 because starts at 0

    save(strcat(output_folder,"/output_table.mat"),"output_table");

    save(strcat(output_folder,"/Rsquared_for_different_factor_nb.mat"),"original_rsquared");
    Rsquared_table = num2cell(original_rsquared);
    rsquared_var_names = 0:numel(original_rsquared)-1;
    Rsquared_table = cell2table(Rsquared_table,"VariableNames",string(rsquared_var_names));

    Rsquared_surrogate_table = num2cell(surrogate_rsquared_mean);
    Rsquared_surrogate_table(2,:) = num2cell(upper_bound); %95th prctile
    Rsquared_surrogate_table(3,:) = num2cell(lower_bound); %5th prctile
    rsquared_var_names = 0:numel(surrogate_rsquared_mean)-1;
    Rsquared_surrogate_table = cell2table(Rsquared_surrogate_table,"VariableNames",string(rsquared_var_names), ...
        "RowNames",{'mean','95th_prctl','5th_prctl'});
    save(strcat(output_folder,"/Surrogate_Rsquared.mat"),"Rsquared_surrogate_table");

    writetable(output_table,strcat(output_folder,"/output_table.csv"));
    writetable(Rsquared_table,strcat(output_folder,"/Rsquared_for_different_factor_nb.csv"));
    writetable(Rsquared_surrogate_table,strcat(output_folder,"/Surrogate_Rsquared.csv"),"WriteRowNames",true);
    
    
end % end of "for each condition"
