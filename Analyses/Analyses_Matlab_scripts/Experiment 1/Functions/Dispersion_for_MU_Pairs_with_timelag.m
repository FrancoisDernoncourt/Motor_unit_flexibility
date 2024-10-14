function [max_dispersion,norm_of_max_dispersion,lag_to_minimize,dispersions,norms] = ...
Dispersion_for_MU_Pairs_with_timelag( ...
mu_one, mu_two, sampling_rate, norm_epsilon, max_timelag, downsampling_factor, plotting, figure_idx)
%1st argument = filtered signal of the first MU of the pair (filtered signal = smoothed discharge rate)
%2nd argument = 2nd MU of the pair
%3rd argument = sampling rate of the signal
%norm epsilon = increments for the norms to be checked

%%

% %For test purpose:
% %     % For condition sin 0.25
% % mu_one = edition.filteredMUs{4}(9,:);
% % mu_two = edition.filteredMUs{6}(8,:);
%     % For condition sin 3
% mu_one = edition.filteredMUs{4}(16,:);
% mu_two = edition.filteredMUs{6}(5,:);
%     %
% sampling_rate = 2048; % in samples per s
% norm_epsilon = 1;
% max_timelag = 0; %in s


% Removing valus which are not smooth (big change of values from one data
% point to the next)
mu_one_diff = abs(diff(mu_one));
mu_two_diff = abs(diff(mu_two));
for i=1:2
    if i == 1
        [~,datapoints_removal] = find(mu_one_diff > 1); % a more than 1hz-gap between data points
    elseif i == 2
        [~,datapoints_removal] = find(mu_two_diff > 1); % a more than 1hz-gap between data points
    end
    if ~isempty(datapoints_removal)
        for j=1:length(datapoints_removal)
            if i == 1
                mu_one(datapoints_removal+1) = NaN;
            elseif i == 2
                mu_two(datapoints_removal+1) = NaN;
            end
        end
    end
end
% mu_one(mu_one==0) = NaN;
% mu_two(mu_two==0) = NaN;

timelag_in_samples = round(max_timelag*sampling_rate);
mu_two_lag = zeros(timelag_in_samples*2 +1 , length(mu_two)-(2*timelag_in_samples) );
mu_one = mu_one(timelag_in_samples+1:end-timelag_in_samples);
for lagi=1:size(mu_two_lag,1)
    lag = lagi-timelag_in_samples-1;
    min_lag = timelag_in_samples + lag + 1;
    max_lag = length(mu_two) - (2*timelag_in_samples + 1) + min_lag;
    temp_mu_two = mu_two(min_lag:max_lag);
    mu_two_lag(lagi,:) = temp_mu_two;
end

%surf(mu_one,1:205,mu_two_lag,'MeshStyle','none');

% Check intersection at norm epsilon for (lag = timelag_in_samples+1) and then, for each intersection, minimize the dispersion over the allowed lags
% Very long, unless timelag is small

MU_pair_norm_total = mu_one + mu_two_lag(timelag_in_samples+1,:);

epsilon_min = max([ceil(min(MU_pair_norm_total)),1]);
epsilon_max = floor(max(MU_pair_norm_total));
% if isnan(epsilon_max)
%     return  
% end

norms = [];
dispersions = cell(length(epsilon_min:1:epsilon_max),1);

if plotting
    if ishandle(figure_idx)
        close(figure(figure_idx));
    end
    figure(figure_idx);
    lim_plot = [0,(max(mu_one)+max(mu_two)+1),...
        0,(max(mu_one)+max(mu_two)+1)]; %1 is x min ; 2 is x max ; 3 is y min ; 4 is y max
%     lim_plot = [0,25,...
%         0,25]; %1 is x min ; 2 is x max ; 3 is y min ; 4 is y max
    % colors = winter((epsilon_max-epsilon_min)+1);
    % CUSTOM COLORMAP
        color1 = [0.2, 0.4, 0.8]; % RGB values for the first color
        color2 = [0.8, 0.2, 0.4]; % RGB values for the second color
        numColors = (epsilon_max-epsilon_min)+1;
        interpValues = linspace(0, 1, numColors);
        customColormap = interp1([0, 1], [color1; color2], interpValues);
        colors = customColormap(1:(epsilon_max-epsilon_min)+1,:);
    plot(mu_one,mu_two_lag(timelag_in_samples+1,:),'linewidth',3,'Color',[0.5,0.5,0.5]);
    title("State-space smoothed discharge rates of MU pair - Dispersion calculation")
    hold on
    scatter(mu_one,mu_two_lag(timelag_in_samples+1,:),3.5,'o','MarkerFaceColor',[0,0,0],'MarkerEdgeAlpha',0);
    axis square
    xlim([lim_plot(1),lim_plot(2)]);
    ylim([lim_plot(3),lim_plot(4)]);
    xlabel('MU X discharge rate (spikes per second)');
    ylabel('MU Y discharge rate (spikes per second)');
end

max_dispersion_for_display = 0;

for i=epsilon_min:1:epsilon_max
    norms(end+1) = i;
    % find all the time points at which the trajectory of mu1 and mu2 in
    % the state-space cross the line defined by norm = i
    dispersions{i} = zeros(1,size(mu_two_lag,1));

    for lagi=1:size(mu_two_lag,1)
        % find the timelag for the max intersec_timepoints distance between 1 and 2 is minimized, and take the max value for
        % this timelag
        MU_pair_norm = mu_one + mu_two_lag(lagi,:);
        % epsilon_range_for_crossing = %in discharge rates = will accept all the intersec_timepoints for which ((MU_pair_norm > i-range) && (MU_pair_norm < i+range))
        if downsampling_factor == 0
            epsilon_range_for_crossing = 0.2;
        else
            epsilon_range_for_crossing = min([0.2 * log(downsampling_factor), norm_epsilon/2]);
        end
        % arbitrary value which seems to work most of the time...

        considered_intersec_timepoints = [];
        temp_epsilon_range_for_crossing = epsilon_range_for_crossing; % gradually increase the acceptable range to prevent error
        while isempty(considered_intersec_timepoints)
            considered_intersec_timepoints = (MU_pair_norm>(i-temp_epsilon_range_for_crossing));
            considered_intersec_timepoints = considered_intersec_timepoints - (MU_pair_norm > (i+temp_epsilon_range_for_crossing) );
            considered_intersec_timepoints = find(considered_intersec_timepoints);
            if ~isempty(considered_intersec_timepoints)
                break;
            else
                temp_epsilon_range_for_crossing = temp_epsilon_range_for_crossing + epsilon_range_for_crossing;
            end
        end
%         if temp_epsilon_range_for_crossing > epsilon_range_for_crossing % just to know the range needed for some epsilons
%             disp(strcat("epislon range for epsilon = ", num2str(i), ...
%                 " at lag = ", num2str(lagi), " : ", num2str(temp_epsilon_range_for_crossing)));
%         end
        indices = diff(considered_intersec_timepoints) ~= 1;  % Find indices where the consecutive sequence breaks
        indices = [true, indices];  % Include the first element
        chunks = mat2cell(considered_intersec_timepoints, 1, diff(find([indices, true]))); % Split the vector into chunks
        intersec_timepoints = cellfun(@median,chunks);
        intersec_timepoints = round(intersec_timepoints);

        % find intersec_timepoints for which the distance between mu1 and 2 is the greatest
        interesc_coord = nan(numel(intersec_timepoints),2);
        for intersections=1:numel(intersec_timepoints)
            interesc_coord(intersections,1) = mu_one(intersec_timepoints(intersections));
            interesc_coord(intersections,2) = mu_two_lag(lagi,intersec_timepoints(intersections));
        end
        [~,first_intersec] = min(interesc_coord(:,1));
        [~,last_intersec] = max(interesc_coord(:,1));
        first_point = interesc_coord(first_intersec,:);
        last_point = interesc_coord(last_intersec,:);

        dispersions{i}(lagi) = sqrt( (last_point(1)-first_point(1))^2 + (last_point(2)-first_point(2))^2 );

        if lagi == timelag_in_samples+1

            % save values for max disp in order to be able to plot it later
            % (only for timelag = +103)
            if plotting % this will plot only for timelag = +102
                lineObj = line([i,0],[0,i],'linestyle','--','color',[colors(i-(epsilon_min-1),:),0.3],'linew',2);
                %lineObj.Color(4) = 0.2;
                %scatter(mu_one(intersec_timepoints),mu_two(intersec_timepoints),70,colors(i-(epsilon_min-1),:),'o','filled');
                line([first_point(1),last_point(1)],[first_point(2),last_point(2)],'color',[colors(i-(epsilon_min-1),:),0.65],'linew',2);
                scatter([first_point(1),last_point(1)],[first_point(2),last_point(2)],200,'x','MarkerEdgeColor',colors(i-(epsilon_min-1),:),'LineWidth',2.5,"MarkerEdgeAlpha",1);
            end
            
            if plotting
                % Just for display with the right time lag
                display_intersec_timepoints = intersec_timepoints;
                display_first_point = first_point;
                display_last_point = last_point;
            end

        end % end of "if no lag (lagi == timelag_in_samples+1)"

    end % end of "lagi"

    current_minimized_dispersion = min(dispersions{i});

    if current_minimized_dispersion > max_dispersion_for_display && plotting
        max_dispersion_for_display = current_minimized_dispersion;
        max_disp_intersec_timepoints = display_intersec_timepoints;
        max_disp_first_points = display_first_point;
        max_disp_last_points = display_last_point;
    end

end % end of "for each norm epsilon"

[temp_max_disp,temp_lag_to_minimize] = cellfun(@min,dispersions,'UniformOutput',false); %minimize across lags, but maximize across dispersion 
temp_max_disp = [temp_max_disp{:}];
temp_lag_to_minimize = [temp_lag_to_minimize{:}] - lag -1;
[max_dispersion,max_dispersion_index] =  max(temp_max_disp);
lag_to_minimize = temp_lag_to_minimize(max_dispersion_index);
norm_of_max_dispersion = norms(max_dispersion_index);
%intersec_timepoints_of_max_dispersion = [1000,5000];

if plotting
    lineObj = line([norm_of_max_dispersion,0],[0,norm_of_max_dispersion],'color','red','linew',3);
    lineObj.Color(4) = 0.3;
    line([max_disp_first_points(1),max_disp_last_points(1)],...
        [max_disp_first_points(2),max_disp_last_points(2)],'color','red','linew',5);
    scatter([max_disp_first_points(1),max_disp_last_points(1)],...
        [max_disp_first_points(2),max_disp_last_points(2)],200,[1 0 0],'o','filled');
end

end % end of function