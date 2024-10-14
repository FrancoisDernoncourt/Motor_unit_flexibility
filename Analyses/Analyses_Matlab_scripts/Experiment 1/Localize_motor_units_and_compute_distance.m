close all
clear all
clc

% load files (not windowed) filtered in discharge rates
Subject_idx = 1;
output_name = strcat("S",num2str(Subject_idx),"_VL_");
condition_names = {"plateau";"sin0.25";"sin1";"sin3"};
corresponding_columns_in_MUs_matched = [4,5,6,7];
nb_of_conditions = numel(corresponding_columns_in_MUs_matched);
grids_to_load = 3:6; %1:4; %1:6
nb_grids = numel(grids_to_load);

% PARAMETERS
% change following parameters according to the muscle (VL or GM)
inter_electrode_distance = 8; %in mm
inter_grid_distance_leftright = 3; % in inter-electrode distance. More than up-down to account for -1 electrode because of differential EMG
inter_grid_distance_updown = 2; % in inter-electrode distance
grids_loc = { "distal medial electrode grid";...
    "distal lateral electrode grid";...
    "proximal medial electrode grid";...
    "proximal lateral electrode grid"};
empty_electrode = [5,12]; % index of the missing electrode channel on each grid
%(take into account that at least +1 electrode because of differential EMG)
    % flip_grid_updown = true;
Electrode_to_discard_from_analysis = [];
fsamp = 2048;

new_dir_name = "MU_pairs_MUAP_distance";
draw_and_save_figures = false;

mkdir(new_dir_name);
mkdir(strcat(new_dir_name,"\plot_output"));
vars_to_save = {"epicenter" ; "epicenter_mm" ; "inter_grid_distance_leftright" ; "inter_grid_distance_updown" ;...
    "MU_pairs_distance_matrix" ; "spike_trigger_average_all_conditions" ; "spike_trigger_average_all_conditions_peak_to_peak" ;...
    "all_grids_concat_conditions" ; "all_grids_concat_average_over_conditions" };

%% LOAD FILES

paths_string_to_load = {};
files_string_to_load = {};

% signal files
for filei=1:nb_of_conditions
    [files_string_to_load{end+1}, paths_string_to_load{end+1}] = ...
        uigetfile(".mat",strcat("Please load the signal file corresponding to condition #",num2str(filei)), ...
        strcat("Please load the signal file corresponding to condition #",num2str(filei)), ...
        "Multiselect","off");
end
% matched file
[files_string_to_load{end+1}, paths_string_to_load{end+1}] = ...
    uigetfile(".mat",strcat("Please load the matched MUs file"), ...
    strcat("Please load the matched MUs file"), ...
    "Multiselect","off");

load(strcat(paths_string_to_load{end},files_string_to_load{end})); % Load MUs matched file

for filei=1:numel(files_string_to_load)
    if filei <= nb_of_conditions % signal files
        disp(strcat("Loading file #",num2str(filei),"/",num2str(numel(files_string_to_load)-1)))
        load(strcat(paths_string_to_load{filei},files_string_to_load{filei}));
        % Remove empty grids from grids to load
        grids_to_load_for_file = grids_to_load;
        for gridi=1:nb_grids
            grididx = grids_to_load(gridi);
            if isempty(edition.Pulsetrain{grids_to_load(gridi)})
                grids_to_load_for_file(find(grids_to_load_for_file==grididx)) = [];
            end
        end
        files_each_condition{filei}.parameters = parameters;
        files_each_condition{filei}.signal = struct();
        files_each_condition{filei}.signal.data = signal.data;
        files_each_condition{filei}.signal.coordinates = signal.coordinates;
        files_each_condition{filei}.signal.EMGmask = signal.EMGmask;
        files_each_condition{filei}.discharge_times = edition.Dischargetimes(grids_to_load_for_file,:);
        files_each_condition{filei}.loaded_grids_idx = grids_to_load_for_file;
        % get correspondance between linear indexing of MUs and indexing
        % based on grid+mu nb
        files_each_condition{filei}.mus_idx = [];
        mu_iter = 0;
        for gridi=1:size(edition.Dischargetimes,1)
            for mui=1:numel(edition.Dischargetimes(gridi,:))
                if isempty(edition.Dischargetimes{gridi,mui})
                    continue
                end
                mu_iter = mu_iter+1;
                files_each_condition{filei}.mus_idx(mu_iter,1:2) = nan;
                files_each_condition{filei}.mus_idx(mu_iter,1) = mu_iter;
                files_each_condition{filei}.mus_idx(mu_iter,2) = gridi*100+mui;
            end
        end
    else % MU match file
        break;
    end
end

clearvars signal parameters edition filei winowi savename grids_to_load_for_file

%% Compute MUAP shapes with spike-trigger average

% 1st dim are MUs ; 2nd dim are conditions ; then there are sub-cells for
% each cell
% average window (1) or peak-to-peak value (2)
spike_trigger_average_window_size = round(0.05*fsamp);
min_amplitude_in_files = []; %actually not necessary, because peak to peak amplitude will always be plotted from 0
max_amplitude_in_files = [];

for conditioni = 1:nb_of_conditions
    disp(strcat("Processing condition #",num2str(conditioni)," / ", num2str(nb_of_conditions)));
    min_amplitude_in_files(conditioni) = inf;
    max_amplitude_in_files(conditioni) = 0;

    corresponding_matched_MUs_file = cell(nb_grids,1);
    parameters_file = files_each_condition{conditioni}.parameters;
    signal_file = files_each_condition{conditioni}.signal;
    discharge_times_file = files_each_condition{conditioni}.discharge_times;
    nb_grids = numel(files_each_condition{conditioni}.loaded_grids_idx);
    SIGNAL = cell(nb_grids,1); % reset signal variable
    spike_trigger_average_file_peak_to_peak = []; % reset spike trigger average

    for gridi = 1:nb_grids
        grididx = files_each_condition{conditioni}.loaded_grids_idx(gridi);
        for channeli = 1:parameters_file.nbelectrodes
            SIGNAL{gridi}{signal_file.coordinates{grididx}(channeli,1), signal_file.coordinates{grididx}(channeli,2)} = ...
                signal_file.data((grididx-1)*parameters_file.nbelectrodes+channeli,:);
            discardChannelsVec(signal_file.coordinates{grididx}(channeli,1), signal_file.coordinates{grididx}(channeli,2)) = ...
                signal_file.EMGmask{grididx}(channeli);
        end

        % Discard channels BEFORE doing differential, and for the
        % differential, check with next to next MU instead of just the
        % next one (if NaN)
        % Discard channels to be discarded (set to nan)
        for ch=1:numel(discardChannelsVec)
            if discardChannelsVec(ch)
                SIGNAL{gridi}{ch} = NaN(1,numel(SIGNAL{gridi}{ch}));
            end
        end

        % get differential signal
        for coli=1:size(SIGNAL{gridi},2)
            temp_diff_col = cell2mat(SIGNAL{gridi}(:,coli));
            temp_diff_col = diff(temp_diff_col);
            SIGNAL_diff_temp(1:size(temp_diff_col,1),coli) = mat2cell(temp_diff_col, ones(size(temp_diff_col,1), 1), size(temp_diff_col,2));
        end

        SIGNAL{gridi} = [];
        SIGNAL{gridi} = SIGNAL_diff_temp;
        SIGNAL{gridi} = rot90(SIGNAL{gridi});
    end
    clearvars SIGNAL_diff_temp discardChannelsVec ch channeli coli temp_diff_col

    % spike_trigger_average_file variable =
    % 1st level are grids (MU per grids)
    % 2nd level are MUs
    % 3rd level are grids (spike trigger average of elecrodes per grid)
    % 4th level are electrodes in the grid
    % 5th level are actual spike trigger average samples
    spike_trigger_average_file = cell(nb_grids,1); %1st level

    for gridi=1:nb_grids %find each spike of each MU for all grids

        grididx = files_each_condition{conditioni}.loaded_grids_idx(gridi);
        idx_of_nonempty_MUs_in_grid_temp = find(cellfun(@(x) ~isempty(x), discharge_times_file(gridi,:)));
        nb_mu_temp = numel(idx_of_nonempty_MUs_in_grid_temp);

        if nb_mu_temp >= 1 % if at least 1 MU in the grid
            % setting the right nb of MUs in the first level of the
            % spike_trigger_average_file variable
            spike_trigger_average_file{gridi} = cell(nb_mu_temp,1);

            for mu_iter=1:nb_mu_temp

                mui = idx_of_nonempty_MUs_in_grid_temp(mu_iter);

                % correspondance
                corresponding_matched_MUs_file{gridi}(mu_iter,1:3) = nan;
                % MU idx (grid*100 + mui)
                corresponding_matched_MUs_file{gridi}(mu_iter,1) = ...
                    files_each_condition{conditioni}.loaded_grids_idx(gridi)*100 + mui;
                % linear index of MU for the current file
                temp_mu_idx_find = find(...
                   files_each_condition{conditioni}.mus_idx(:,2)==corresponding_matched_MUs_file{gridi}(mu_iter,1));
                corresponding_matched_MUs_file{gridi}(mu_iter,2) = temp_mu_idx_find;
                % find corresponding MU in 'MUs matched'
                temp_mu_idx_find = find(...
                    [MUs_matched{:,corresponding_columns_in_MUs_matched(conditioni)}{:}]==corresponding_matched_MUs_file{gridi}(mu_iter,1));
                if ~isempty(temp_mu_idx_find)
                    if numel(temp_mu_idx_find) > 1 % in case the sme MU has been matched twice (happened on rare occasions)
                        corresponding_matched_MUs_file{gridi}(mu_iter,3) = temp_mu_idx_find(1); % arbitrarily choose the 1st MU
                    else
                        corresponding_matched_MUs_file{gridi}(mu_iter,3) = temp_mu_idx_find;
                    end
                end
                % Turn into a table if last iteration through MUs
                if mu_iter==nb_mu_temp
                    corresponding_matched_MUs_file{gridi} = array2table(corresponding_matched_MUs_file{gridi},"VariableNames",...
                        {'MU_idx','linear_idx','matched_idx'});
                    clearvars temp_mu_idx_find
                end

                % setting all the sub-levels for the
                % spike_trigge_average_file variable
                spike_trigger_average_file{gridi}{mui} = cell(4,1);
                for gridi2 = 1:nb_grids
                    spike_trigger_average_file{gridi}{mui}{gridi2} = cell(5,12);
                    for electroderowi=1:5
                        for electrodecoli=1:12 % 12 because diff() removed one column
                            spike_trigger_average_file{gridi}{mui}{gridi2}{electroderowi,electrodecoli} = ...
                                zeros(1,(spike_trigger_average_window_size*2)+1);
                        end
                    end
                end

                spike_times = discharge_times_file{gridi,mu_iter};
                if ~isempty(spike_times)

                    for gridi2 = 1:nb_grids
                        for channelrowi = 1:5
                            for channelcoli = 1:12 % 12 because diff() removed one column
                                if isempty(SIGNAL{gridi2}{channelrowi,channelcoli})
                                    spike_trigger_average_file{gridi}{mui}{gridi2}{channelrowi,channelcoli} = nan;
                                    spike_trigger_average_file_peak_to_peak{gridi}{mui}{gridi2}(channelrowi,channelcoli) = nan;
                                    continue;
                                else

                                    for spikei=1:numel(spike_times)
                                        spike_time = spike_times(spikei);
                                        spike_trigger_signal(spikei,:) = SIGNAL{gridi2}{channelrowi,channelcoli}...
                                            (spike_time-spike_trigger_average_window_size:spike_time+spike_trigger_average_window_size);
                                    end % end of "for each spike"

                                    spike_trigger_average_file{gridi}{mui}{gridi2}{channelrowi,channelcoli} = mean(spike_trigger_signal,1,'omitnan');

                                    spike_trigger_average_file_peak_to_peak{gridi}{mui}{gridi2}(channelrowi,channelcoli) = abs( ...
                                        min(spike_trigger_average_file{gridi}{mui}{gridi2}{channelrowi,channelcoli})) + abs(max( ...
                                        spike_trigger_average_file{gridi}{mui}{gridi2}{channelrowi,channelcoli}));

                                end
                            end % end of "for each column of channel / electrode"
                        end % end of "for each row of channel / electrode"

                        % Interpolate missing values (discarded channels)
                        interpolated_mat_temp = spike_trigger_average_file_peak_to_peak{gridi}{mui}{gridi2};
                        [X, Y] = meshgrid(1:size(interpolated_mat_temp,2), 1:size(interpolated_mat_temp,1));
                        nonNanIndices_temp = ~isnan(interpolated_mat_temp);
                        interpolated_mat_temp(~nonNanIndices_temp) = griddata(...
                            X(nonNanIndices_temp), Y(nonNanIndices_temp),...
                            interpolated_mat_temp(nonNanIndices_temp),...
                            X(~nonNanIndices_temp), Y(~nonNanIndices_temp));
                        spike_trigger_average_file_peak_to_peak{gridi}{mui}{gridi2} = interpolated_mat_temp;
                        
                        % set min and max amplitude within the file (for
                        % appropriate scaling in the figures)
                        if max(spike_trigger_average_file_peak_to_peak{gridi}{mui}{gridi2}(:)) > max_amplitude_in_files(conditioni)
                            max_amplitude_in_files(conditioni) = max(spike_trigger_average_file_peak_to_peak{gridi}{mui}{gridi2}(:));
                        end
                        if min(spike_trigger_average_file_peak_to_peak{gridi}{mui}{gridi2}(:)) < min_amplitude_in_files(conditioni)
                            min_amplitude_in_files(conditioni) = min(spike_trigger_average_file_peak_to_peak{gridi}{mui}{gridi2}(:));
                        end

                    end % end of "for each grid" (electrode per grid)
                end % end of "if the MU is firing" (spike count > 0)
            end % end of mui (each matched mu)
        end % end of "if nb of matched MUs in the grid > 1"
    end % end of gridi (MU per grid)
    files_each_condition{conditioni}.spike_trigger_average_peak_to_peak_amplitude = spike_trigger_average_file_peak_to_peak;
    files_each_condition{conditioni}.min_p2p_amplitude = min_amplitude_in_files(conditioni);
    files_each_condition{conditioni}.max_p2p_amplitude = max_amplitude_in_files(conditioni);
    files_each_condition{conditioni}.MU_correspondance_table = corresponding_matched_MUs_file;
end % end of "for each condition / file"

clearvars spike_trigger_signal spike_trigger_average_file spike_trigger_average_file_peak_to_peak...
    nonNanIndices_temp interpolated_mat_temp corresponding_matched_MUs_file

%% GET MUAP AMPLITUDE AVERAGE AND EPICENTER FOR EACH MU

% maybe interesting to convolve average values with a smoothing kernel ?
% Define the convolution kernel (blurring effect)
kernel = [1 1 1; 1 1 1; 1 1 1]./9; % This is a simple averaging kernel

close all
for conditioni=1:nb_of_conditions
    mu_iter = 0;
    files_each_condition{conditioni}.smoothed_MUAP_amplitude = {}; % 1 cell per MU
    % 4 sub-cells per cell, with each sub-cell containing a 5*12 matrix
    files_each_condition{conditioni}.epicenter = []; % 1 cell per MU
    % epicenter = grid#, x pos# (rows of electrodes), y pos# (columns of electrodes)
    for gridi=1:numel(files_each_condition{conditioni}.loaded_grids_idx)

        nb_of_MUs = size(files_each_condition{conditioni}.spike_trigger_average_peak_to_peak_amplitude{gridi},2);
        for mui = 1:nb_of_MUs
            mu_iter = mu_iter+1;
            current_peak_to_peak_amplitude_map = files_each_condition{conditioni}.spike_trigger_average_peak_to_peak_amplitude{gridi}{mui};
            nb_grids_for_epicenter = numel(current_peak_to_peak_amplitude_map);

            epicenter_each_grid_temp = cell(1,nb_grids_for_epicenter);
            max_amplitude_each_grid_temp = zeros(1,nb_grids_for_epicenter);
            true_epicenter_temp = [nan,nan,nan]; %grid#, x pos# (rows of electrodes), y pos# (columns of electrodes)

            temp_smoothed_grid = {};
            if draw_and_save_figures
                figure(conditioni*100+mu_iter)
            end
            for gridi2=1:nb_grids_for_epicenter
                % temp_smoothed_grid{gridi2} = current_peak_to_peak_amplitude_map{gridi2};
                temp_smoothed_grid{gridi2} = conv2(current_peak_to_peak_amplitude_map{gridi2},kernel,'same');
                % Restore missing electrode in the corner
                temp_smoothed_grid{gridi2}(empty_electrode(1),empty_electrode(2)) = nan;
                % Find max amplitude within this grid
                [maxAmplitude_temp, linearIndexOfMaxAmplitude_temp] = max([temp_smoothed_grid{gridi2}(:)]);
                [epicenter_each_grid_temp{gridi2}(1),epicenter_each_grid_temp{gridi2}(2)] = ...
                    ind2sub(size(temp_smoothed_grid{gridi2}),linearIndexOfMaxAmplitude_temp);
                % Check if this is the grid with the current highest
                % amplitude amongst grids
                max_amplitude_each_grid_temp(gridi2) = maxAmplitude_temp;

                % PLOTTING
                if draw_and_save_figures
                    subplot(2,2,gridi2)
                    % define grid
                    [X, Y] = meshgrid(1:size(temp_smoothed_grid{gridi2},2), 1:size(temp_smoothed_grid{gridi2},1));
                    surf_plot = surf(X,Y,temp_smoothed_grid{gridi2});
                    % surf_plot.EdgeColor = 'none';
                    surf_plot.EdgeColor = 'interp';
                    surf_plot.LineWidth = 1; 
                    surf_plot.FaceAlpha = 0.5;
                    hold on
                    for xi = 1:size(temp_smoothed_grid{gridi2},2)
                        for yi = 1:size(temp_smoothed_grid{gridi2},1)
                            scatter3(xi,yi,temp_smoothed_grid{gridi2}(yi,xi),30,...
                                temp_smoothed_grid{gridi2}(yi,xi),'filled');
                        end
                    end
                    % surf_plot.FaceColor = 'interp';
                    caxis([0,max_amplitude_in_files(conditioni)]);
                    zlim([0,max_amplitude_in_files(conditioni)]);
                    xlim([0,13]);
                    ylim([0,6]);
                    colormap turbo;
                    title(strcat("Grid#",num2str(gridi2)," - ",grids_loc{gridi2}))
                    box off
                    ax = gca;
                    ax.XTick = 1:12;
                    ax.YTick = 1:5;
                    % ax.ZTick = [];
                    % ax.XColor = 'none';
                    % ax.YColor = 'none';
                    ax.DataAspectRatio = [1, 1, max_amplitude_in_files(conditioni)/3]; % Set the data aspect ratio to be equal
                    zlabel("Amplitude (\muV)")
                    xlabel("Channel columns")
                    ylabel("Channel rows")
                end
            end

            [maxAmplitude_temp,true_epicenter_temp(1)] = max(max_amplitude_each_grid_temp);
            % reverse X and Y values
            true_epicenter_temp(3) = epicenter_each_grid_temp{true_epicenter_temp(1)}(1);
            true_epicenter_temp(2) = epicenter_each_grid_temp{true_epicenter_temp(1)}(2);

            files_each_condition{conditioni}.smoothed_MUAP_amplitude{mu_iter} = temp_smoothed_grid;
            files_each_condition{conditioni}.epicenter(mu_iter,:) = true_epicenter_temp;

            if draw_and_save_figures
                sgtitle({"MUAP peak to peak amplitude at each channel";...
                    strcat("MU#",num2str(...
                    files_each_condition{conditioni}.MU_correspondance_table{gridi}{mui,1}), ...
                    " (linear index = ",num2str(...
                    files_each_condition{conditioni}.MU_correspondance_table{gridi}{mui,2}),...
                    " ; matched idx = ",num2str(...
                    files_each_condition{conditioni}.MU_correspondance_table{gridi}{mui,3}),")")});
                subplot(2,2,true_epicenter_temp(1))
                hold on
                text(true_epicenter_temp(2)+0.5, true_epicenter_temp(3)-0.25, max_amplitude_in_files(conditioni),...
                    strcat("Grid #", num2str(true_epicenter_temp(1)), "; row:",num2str(true_epicenter_temp(3)),...
                    " - col:", num2str(true_epicenter_temp(2)) ), ...
                    'FontSize', 12, 'FontWeight', 'bold' , 'Color', 'red', ...
                    'BackgroundColor','white');
                line([true_epicenter_temp(2),true_epicenter_temp(2)],...
                    [true_epicenter_temp(3),true_epicenter_temp(3)],...
                    [0,max_amplitude_in_files(conditioni)],...
                    'Color','red','LineWidth',3,'LineStyle',':');
                scatter3([true_epicenter_temp(2),true_epicenter_temp(2)],...
                    [true_epicenter_temp(3),true_epicenter_temp(3)],...
                    [max_amplitude_in_files(conditioni),max_amplitude_in_files(conditioni)],...
                    100,"MarkerFaceColor","white","MarkerEdgeColor","red",...
                    "LineWidth",3);
                scatter3([true_epicenter_temp(2),true_epicenter_temp(2)],...
                    [true_epicenter_temp(3),true_epicenter_temp(3)],...
                    [maxAmplitude_temp,maxAmplitude_temp],...
                    500,"x","MarkerFaceColor","white","MarkerEdgeColor","red",...
                    "LineWidth",5);
                saveas(gcf(),strcat(new_dir_name,"\plot_output\Condition_",num2str(conditioni),"_MU_",num2str(mu_iter),".png"))
                close(gcf())
            end

        end % end of "for each MU" (mui)
    end % end of "for each grid' (gridi)
end % end of "for each file/condition" (conditioni)

clearvars temp_smoothed_grid epicenter_each_grid_temp maxAmplitude_temp nb_mu_temp...
current_peak_to_peak_amplitude_map true_epicenter_temp maxAmplitude_temp linearIndexOfMaxAmplitude_temp

%% GET LOCALIZATION OF EACH MU & DISTANCE BETWEEN PAIRS OF MUS

MU_distances_output_table = [];
distance_table_varnames = {'Subject','Condition','MU_X_idx','MU_X_idx_matched','MU_Y_idx','MU_Y_idx_matched','Distance'};
dist_list = [];

for conditioni=1:nb_of_conditions
    files_each_condition{conditioni}.output_table = files_each_condition{conditioni}.MU_correspondance_table{1};
    for gridi=2:numel(files_each_condition{conditioni}.MU_correspondance_table)
        files_each_condition{conditioni}.output_table = vertcat(...
            files_each_condition{conditioni}.output_table, ...
            files_each_condition{conditioni}.MU_correspondance_table{gridi});
    end
    MU_list_position = array2table(files_each_condition{conditioni}.epicenter,"VariableNames",...
        {'Grid','Column_channel','Row_channel'});
    files_each_condition{conditioni}.output_table = horzcat(...
        files_each_condition{conditioni}.output_table, MU_list_position);
    
    % Get position in mm
    % Starting from grid 3, channel at row 1 & column 1 (bottom left of the 4 grids)
    MU_list_position = [];
    nb_MUs = size(files_each_condition{conditioni}.output_table,1);
    for mui=1:nb_MUs
        MU_list_position(mui,1:2) = nan;
        grid_temp = files_each_condition{conditioni}.output_table.Grid(mui);

        go_right = false;
        go_up = false;
        offset_row = 0;
        offset_col = 0;
        switch grids_loc{grid_temp}
            case "proximal medial electrode grid" % UP from origin (grid 1 for VL)
                go_up = true;
            case "distal medial electrode grid" % UP and RIGHT from origin (grid 2 for VL)
                go_up = true;
                go_right = true;
            case "proximal lateral electrode grid" % ORIGIN (grid 3 for VL)
            case "distal lateral electrode grid" % RIGHT from origin (grid 4 for VL)
                go_right = true;
        end
        if go_up
            offset_row = 5 + inter_grid_distance_updown;
        end
        if go_right
            offset_col = 12 + inter_grid_distance_leftright;
        end
        MU_list_position(mui,1) = (files_each_condition{conditioni}.output_table.Row_channel(mui) + offset_row)*inter_electrode_distance;
        MU_list_position(mui,2) = (files_each_condition{conditioni}.output_table.Column_channel(mui) + offset_col)*inter_electrode_distance;
    end
    MU_list_position = array2table(MU_list_position,"VariableNames",...
        {'latero_medial_pos','proxi_distal_pos'});
    files_each_condition{conditioni}.output_table = horzcat(...
        files_each_condition{conditioni}.output_table, ...
        MU_list_position);
    MU_positions_table = files_each_condition{conditioni}.output_table;
    save(strcat(...
        new_dir_name, "/", output_name, "condition_",num2str(conditioni),"_MU_positions.mat"), ...
        'MU_positions_table');

    mu_pair_iter = size(dist_list,1);
    for mux=1:nb_MUs-1
        for muy=mux+1:nb_MUs
            mu_pair_iter = mu_pair_iter + 1;
            dist_list(mu_pair_iter,:) = nan(1,numel(distance_table_varnames));
            % {'Subject','Condition','MU_X_idx','MU_X_idx_matched','MU_Y_idx','MU_Y_idx_matched','Distance'};
            dist_list(mu_pair_iter,1) = Subject_idx;
            dist_list(mu_pair_iter,2) = conditioni;
            dist_list(mu_pair_iter,3) = files_each_condition{conditioni}.output_table.linear_idx(mux);
            dist_list(mu_pair_iter,4) = files_each_condition{conditioni}.output_table.matched_idx(mux);
            dist_list(mu_pair_iter,5) = files_each_condition{conditioni}.output_table.linear_idx(muy);
            dist_list(mu_pair_iter,6) = files_each_condition{conditioni}.output_table.matched_idx(muy);

            pos_mux = MU_list_position{mux,:};
            pos_muy = MU_list_position{muy,:};
            dist_temp = sqrt(...
                (pos_muy(1)-pos_mux(1))^2 + ...
                (pos_muy(2)-pos_mux(2))^2 );

            dist_list(mu_pair_iter,7) = dist_temp;
        end
    end % end of looping through MU pairs
end % end of "for each condition"
dist_list(isnan(dist_list)) = 0;

MU_distances_output_table = array2table(dist_list,"VariableNames",distance_table_varnames);
save(strcat(...
        new_dir_name, "/", output_name, "distances_between_MUs.mat"), ...
        'MU_distances_output_table');
writetable(MU_distances_output_table,strcat(...
        new_dir_name, "/", output_name, "distances_between_MUs.csv"));

clearvars MU_positions_table mux muy 
