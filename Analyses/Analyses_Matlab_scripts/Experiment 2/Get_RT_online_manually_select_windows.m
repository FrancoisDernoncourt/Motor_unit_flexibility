close all
clear all
clc

%% Get file
[filename,path] = ...
    uigetfile(".mat",strcat("Select the baseline edited file on which to calculate RTs"), ...
    strcat("Select the baseline edited file on which to calculate RTs"), ...
    "Multiselect","off");

[filename_force,path_force] = ...
    uigetfile(".mat",strcat("Select file with force values"), ...
    strcat("Select file with force values"), ...
    "Multiselect","off");


%% Parameters
mean_or_median_RT = "mean"; %"mean
nb_of_trials = 4;
first_spike_burst_max_time = 1; % in seconds - maximum time interval to have the 3 first spikes
fsamp = 2048;
% when this condition is met, the first spike time of the series will be
% considered as the time of recruitment

disp("Loading files");
loaded_file = strcat(path,filename);
load(loaded_file);
loaded_file = strcat(path_force,filename_force);
load(loaded_file);

% Convert force signal into MVC %
% reverse automatically if force reversed = and remove beginning and end to
% avoid artifacts
if abs(min(signal.path(fsamp:end-fsamp))) > abs(max(signal.path(fsamp:end-fsamp)))
    force_percent = signal.path .* -1;
else
    force_percent = signal.path;
end
%force_percent = force_percent - force_baseline;
force_percent = force_percent ./ force_MVC;
force_percent = force_percent .* 100;
% min(force_percent)
% max(force_percent)

%% Select trials
start_of_trials = [];
end_of_trials = [];
for triali = 1:nb_of_trials
    figure()
    plot(force_percent)
    title(strcat("Please select the edges of the window for trial #",num2str(triali)))
    [x_select,~] = ginput(2);
    start_of_trials(triali) = x_select(1);
    end_of_trials(triali) = x_select(2);
    close(figure)
end
start_of_trials = round(start_of_trials);
end_of_trials = round(end_of_trials);

%% Get all trials and the corresponding samples
samples_to_consider_each_trial = cell(1,nb_of_trials);
for triali = 1:nb_of_trials
    samples_to_consider_each_trial{triali} = start_of_trials(triali):end_of_trials(triali);
end

%% Get binary discharge matrix
MUbinary_mat = zeros(1,length(edition.time));
corresponding_MU_in_mat_rows = zeros(1,2);
mu_iter = 0;
for gridi = 1:size(edition.Pulsetrain,2)
    for mui = 1:size(edition.Pulsetrain{gridi},1)
        mu_iter = mu_iter + 1;
        corresponding_MU_in_mat_rows(mu_iter,:) = [gridi,mui];
        spike_times_temp = edition.Dischargetimes{gridi,mui};
        MUbinary_mat(mu_iter,:) = zeros(1,length(edition.time));
        MUbinary_mat(mu_iter,spike_times_temp) = 1;
    end
end

%% GET RECRUITMENT THRESHOLDS
% Rule for Recruitment Time = first spike of the first series of 3 spikes in a row (3
% spikes in < first_spike_burst_max_time)

RT_per_trial = nan(size(MUbinary_mat,1),nb_of_trials);
for triali = 1:nb_of_trials
    binary_mat_temp = MUbinary_mat(:,samples_to_consider_each_trial{triali});
    force_for_trial_temp = force_percent(samples_to_consider_each_trial{triali});
    first_spike_per_mu_for_this_trial = [];
    for mui=1:size(binary_mat_temp,1)
        spike_idx_temp = find(binary_mat_temp(mui,:));
        % if not enough spike, dismiss (less than 4 spikes)
        if numel(spike_idx_temp) < 3
            RT_per_trial(mui,triali) = nan;
            first_spike_per_mu_for_this_trial(mui) = nan;
            continue;
        end
        % take the first spike of the first burst of at least 3
        % spikes happening in < first_spike_burst_max_time
        if isempty(spike_idx_temp)
            RT_per_trial(mui,triali) = nan;
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
            RT_per_trial(mui,triali) = nan;
        else
            RT_per_trial(mui,triali) = force_for_trial_temp(first_spike_per_mu_for_this_trial(mui));
        end
    end % end of "for each MU"

    % Visual check
    close(figure(triali+1))
    figure(triali+1)
    % FORCE SIGNAL
    plot(force_for_trial_temp,'LineWidth',1.5);
    hold on
    plot(force_for_trial_temp,'LineWidth',1.5);
    for mui=1:length(first_spike_per_mu_for_this_trial)
        if mui > length(first_spike_per_mu_for_this_trial) % if the MU doesn't discharge during this trial
            break;
        end
        color_temp = [.5,.5,.5];
        if first_spike_per_mu_for_this_trial(mui) > 0
            scatter(first_spike_per_mu_for_this_trial(mui),mui,'o','MarkerEdgeColor','red');
        end
        spike_train_plot_temp = binary_mat_temp(mui,:);
        spike_train_plot_temp(spike_train_plot_temp==0) = nan;
        spike_train_plot_temp = spike_train_plot_temp .* mui;
        plot(spike_train_plot_temp,'|','Color',color_temp);
    end
    title(strcat("Sanity check for recruitment thresholds - trial #",num2str(triali)));
    saveas(gcf(),strcat("Recruitment - trial ",num2str(triali),".png"));
    % end of visual check

end % end of "for each trial"

switch mean_or_median_RT
    case "mean"
        recruitment_thresholds = mean(RT_per_trial,2,'omitnan'); %in % MVC
    case "median"
        recruitment_thresholds = median(RT_per_trial,2,'omitnan'); %in % MVC
    otherwise
        recruitment_thresholds = mean(RT_per_trial,2,'omitnan'); %in % MVC
end
recruitment_thresholds = [corresponding_MU_in_mat_rows(:,1)*100 + corresponding_MU_in_mat_rows(:,2), recruitment_thresholds];

%%

vars_to_save = {"recruitment_thresholds"};
save(strcat("Recruitment_thresholds_",mean_or_median_RT),vars_to_save{:});

