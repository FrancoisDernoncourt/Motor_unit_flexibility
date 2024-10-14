
%% PARAMETERS

close all
clear all

plot_MUAPs = false;

output_dir_name = "MATCHING_output";
vars_to_save = {"matches_table" ; "MUs_matched" ; "nb_of_MU_each_file" ; "max_nb_of_matched_MUs" };
plot_vars_to_save = {"MUAP_window_size" ; "MUAP_shapes_to_save"};

%% Load files to match

number_of_files_to_match = 4;
for n=1:number_of_files_to_match
    [filename,pathname] = uigetfile('*.mat', 'Select edited file (Avrillon format)', 'MultiSelect', 'off');
    files{n,1} = pathname;
    files{n,2} = filename;
end
clearvars filename pathname n

%% MATCHING !

% for 'Dischargetimes' and 'Pulsetrains' =
% Use data.edition if already processed, otherwise just data.signal

nb_of_MU_each_file = {};

for n = 1:size(files,1)-1
    for m = n+1:size(files,1)

        data1 = load(strcat(files{n,1},files{n,2})); % Load file 1
        data2 = load(strcat(files{m,1},files{m,2})); % Load file 2

        nmu1 = 1;
        nmu2 = 1;
        PulseT1 = [];
        PulseT2 = [];
        Distime1 = {};
        Distime2 = {};

        for j = 1:data1.signal.ngrid % Loop on each grid
            EMGmask = data1.signal.EMGmask{j} + data2.signal.EMGmask{j}; % Get the masks for boths grids (we remove all the channels masked on either task 1 or task 2
            for k = 1:size(data1.edition.Pulsetrain{j},1) % Format the discharge times for the 1st grid
                Distim{k} = data1.edition.Dischargetimes{j,k};
                idxmu1(nmu1) = j*100+k; % Use an index for each motor unit (e.g., 101 is the first motor unit of grid 1)
                nmu1 = nmu1 + 1;
            end
            if ~isempty(k)
                datatmp1 = data1.signal.data((j-1)*64+1:(j-1)*64+length(EMGmask),:);
                datatmp1 = bandpassingals(datatmp1,data1.signal.fsamp, 1);
                datatmp2 = data2.signal.data((j-1)*64+1:(j-1)*64+length(EMGmask),:);
                datatmp2 = bandpassingals(datatmp2,data2.signal.fsamp, 1);
                MUFilters = getMUfilters(datatmp1, EMGmask, Distim); % Get the MU filters from the first file
                [PulseT, Distime, ~] = getPulseT(datatmp2, EMGmask, MUFilters, data1.signal.fsamp);  % Apply the MU filters on the second file
                PulseT1 = [PulseT1; PulseT]; % Concatenate all the grids (file 1)
                Distime1 = [Distime1, Distime];
            end
            clearvars Distime PulseT Distim MUFilters EMGmask datatmp1 datatmp2
        end

        for j = 1:data2.signal.ngrid % Loop on each grid
            PulseT2 = [PulseT2; data2.edition.Pulsetrain{j}]; % Concatenate all the grids (file 2)
            for k = 1:size(data2.edition.Pulsetrain{j},1)
                Distime2{nmu2} = data2.edition.Dischargetimes{j,k};
                idxmu2(nmu2) = j*1000+k; % Use an index for each motor unit (e.g., 1001 is the first motor unit of grid 1)
                nmu2 = nmu2+1;
            end
        end

        nb_of_MU_each_file{1,end+1} = files{n,2};
        MU_per_grid_counting = cellfun(@size, data1.edition.Pulsetrain, UniformOutput=false);
        nb_of_MU_each_file{2,end} = 0;
        for i=1:numel(MU_per_grid_counting)
            nb_of_MU_each_file{2,end} = nb_of_MU_each_file{2,end} + MU_per_grid_counting{i}(1);
        end
        nb_of_MU_each_file{1,end+1} = files{m,2};
        MU_per_grid_counting = cellfun(@size, data2.edition.Pulsetrain, UniformOutput=false);
        nb_of_MU_each_file{2,end} = 0;
        for i=1:numel(MU_per_grid_counting)
            nb_of_MU_each_file{2,end} = nb_of_MU_each_file{2,end} + MU_per_grid_counting{i}(1);
        end
        clearvars MU_per_grid_counting

        Pulsetemp = [PulseT2; PulseT1]; % Concatenate the Pulse trains generated with the filters of the file 1 and the Pulse trains from file 2
        Distemp = [Distime2 Distime1];
        idxmus = [idxmu2 idxmu1];
        [Pulsetemp, Distemp, matchetemp] = check_duplicates(Pulsetemp, Distemp, idxmus, round(data2.signal.fsamp/40), 0.00025, 0.3, data2.signal.fsamp); % Find the matches

        % Here you will need to sort the matches and the new Pulse trains
        % 2 values in a row (matches between files 1 and 2)
        % 1 value in the hundreds (New pulse train without a match), if editable,
        % new MU and new match
        % 1 value in the thousands, no match
        [~, idx] = sort(matchetemp(:,2), 'descend');

        for i = 1:length(idx)
            matches{m, n}(i,:) = matchetemp(idx(i), :);
            PulseT{m, n}(i,:) = Pulsetemp(idx(i), :);
            distimenew{m, n}{i} = Distemp{idx(i)};
        end
        clearvars -except files n m matches PulseT distimenew...
            nb_of_MU_each_file number_of_files_to_match plot_MUAPs...
            output_dir_name vars_to_save plot_vars_to_save
    end
end

[~, unique_indices] = unique({nb_of_MU_each_file{1,:}},'stable');
nb_of_MU_each_file = nb_of_MU_each_file(:,unique_indices);
clearvars unique_indices

%% OUTPUT

% HOW TO INTERPRET
% XXXX = MU from file 2 (file from the ROW)
% XXX = MU from file 1 (file from the COLUMN)
% [XXXX,0] = no match
% [XXX,0] = potential match (new pulse train) if editable
% [XXXX,XXX] = match between file 1 (MU XXXX) and 2 (MU XXX)

%%%%%%%%%%%%% IF NEEDED, REDUCE FILENAME LENGTH HERE
%reduce_filename = numel(['_MUs_selected_and_signal_filtered_DR_Lowpass.mat']);
    % for filei=1:size(files,1)
    %     temp_filename = files(filei,2);
    %     files{filei,2} = temp_filename{:}(1:end-reduce_filename);
    % end

matches_table = cell2table(matches,'VariableNames',files(1:end-1,2),'RowNames',files(:,2));

max_nb_of_MU = max([nb_of_MU_each_file{2,:}]);
MUs_matched = cell(max_nb_of_MU,size(files,1)+4);
MUs_matched = cell2table(MUs_matched,'VariableNames',[{'MU_idx'};{'Muscle'};{'Grid'};files(:,2);'MU_found_in_how_many_files']);

for colmatch=1:size(matches_table,2)
    for rowmatch=1:size(matches_table,1)
        if ~cellfun(@isempty,matches_table{rowmatch,colmatch})
            row_filename = matches_table.Properties.RowNames(rowmatch);
            col_filename = matches_table.Properties.VariableNames(colmatch);

            thousands_file_MUs = matches_table{rowmatch,colmatch}{:}(:,1);
            hundreds_file_MUs = matches_table{rowmatch,colmatch}{:}(:,2);
            %REMOVE "zeros" OR "hundred" ROWS IN
            % thousands files OR rows which have both thousands on either side
            % and then remove the same rows in
            % "hundreds_file"
            idx_to_remove = [find(thousands_file_MUs<1000);find(hundreds_file_MUs>1000)];
            idx_to_remove = sort(idx_to_remove,'descend');
            thousands_file_MUs(idx_to_remove) = [];
            hundreds_file_MUs(idx_to_remove) = [];

            for mu_idx_in_file = 1:numel(MUs_matched{:,row_filename})
                if mu_idx_in_file > numel(thousands_file_MUs)
                    break;
                elseif thousands_file_MUs(mu_idx_in_file) < 1 || hundreds_file_MUs(mu_idx_in_file) < 1
                    continue;
                end
                %%%%%%%%%%%%% If the code goes up to here, there is a match %%%%%%%%%%%%%%%%%

                % grid of identified MU may be different in different files
                % (for VL at least)
                mu_grid_1 = floor(thousands_file_MUs(mu_idx_in_file)/1000);
                mu_grid_2 = floor(hundreds_file_MUs(mu_idx_in_file)/100);
                mu_idx = mod(thousands_file_MUs(mu_idx_in_file),1000);
                mu_curently_considered_1 = mu_grid_1*100+mu_idx;
                mu_curently_considered_2 = hundreds_file_MUs(mu_idx_in_file);

                % Case 1 : match with a new MU to create in the table => go
                % up to the first empty row
                if sum([MUs_matched{:,row_filename}{:}]==mu_curently_considered_1) < 1 ...
                        && sum([MUs_matched{:,col_filename}{:}]==mu_curently_considered_2) < 1
                    % Find the first empty row and fill it
                    first_empty_row = 1;
                    empty_row_found = false;
                    for rowi=1:size(MUs_matched,1)
                        row_content = [MUs_matched{rowi,:}{:}];
                        if isempty(row_content)
                            first_empty_row = rowi;
                            empty_row_found = true;
                            break;
                        else
                            continue;
                        end
                    end
                    if empty_row_found
                        MUs_matched{first_empty_row,row_filename} = {mu_curently_considered_1};
                        MUs_matched{first_empty_row,col_filename} = {mu_curently_considered_2};
                        MUs_matched{first_empty_row,"Grid"}{:}(end+1:end+2) = [mu_grid_1,mu_grid_2];
                    else
                        warning("No empty row found = it means that there is a conflict somewhere");
                    end

                    % Case 2 : MU2 matches with an already existing MU in the table
                elseif sum([MUs_matched{:,row_filename}{:}]==mu_curently_considered_1) < 1 ...
                        && sum([MUs_matched{:,col_filename}{:}]==mu_curently_considered_2) >= 1
                    % FInd the row of the match of MU2
                    similar_idx = 1;
                    for iterating_through_table_rows=1:numel(MUs_matched{:,col_filename})
                        row_content = [MUs_matched{iterating_through_table_rows,col_filename}{:}];
                        if sum(row_content==mu_curently_considered_2) >= 1
                            similar_idx = iterating_through_table_rows;
                            break;
                        else
                            continue;
                        end
                    end
                    MUs_matched{similar_idx,row_filename} = {mu_curently_considered_1};
                    MUs_matched{similar_idx,col_filename}= {[MUs_matched{similar_idx,col_filename}{:},mu_curently_considered_2]};
                    MUs_matched{similar_idx,"Grid"}{:}(end+1:end+2) = [mu_grid_1,mu_grid_2];
                    % Add the MUs into the column of the other file (the
                    % "column file")

                    % Case 3 : MU1 matches with an already existing MU in the table
                elseif sum([MUs_matched{:,row_filename}{:}]==mu_curently_considered_1) >= 1 ...
                        && sum([MUs_matched{:,col_filename}{:}]==mu_curently_considered_2) < 1
                    % FInd the row of the match of MU1
                    similar_idx = 1;
                    for iterating_through_table_rows=1:numel(MUs_matched{:,row_filename})
                        row_content = [MUs_matched{iterating_through_table_rows,row_filename}{:}];
                        if sum(row_content==mu_curently_considered_1) >= 1
                            similar_idx = iterating_through_table_rows;
                            break;
                        else
                            continue;
                        end
                    end
                    MUs_matched{similar_idx,row_filename}= {[MUs_matched{similar_idx,row_filename}{:},mu_curently_considered_1]};
                    MUs_matched{similar_idx,col_filename} = {mu_curently_considered_2};
                    MUs_matched{similar_idx,"Grid"}{:}(end+1:end+2) = [mu_grid_1,mu_grid_2];

                    % Case 4 : MU1 & MU2 match with already existing MUs in the
                    % table
                    % in the case, prioritize MU1 to fill the columns. If there
                    % is a conflict somewhere, it will be seen later.
                elseif sum([MUs_matched{:,row_filename}{:}]==mu_curently_considered_1) >= 1 ...
                        && sum([MUs_matched{:,col_filename}{:}]==mu_curently_considered_2) >= 1
                    % FInd the row of the match of MU1
                    similar_idx_1 = 1;
                    for iterating_through_table_rows=1:numel(MUs_matched{:,row_filename})
                        row_content = [MUs_matched{iterating_through_table_rows,row_filename}{:}];
                        if sum(row_content==mu_curently_considered_1) >= 1
                            similar_idx_1 = iterating_through_table_rows;
                            break;
                        else
                            continue;
                        end
                    end
                    % FInd the row of the match of MU2
                    similar_idx_2 = 1;
                    for iterating_through_table_rows=1:numel(MUs_matched{:,col_filename})
                        row_content = [MUs_matched{iterating_through_table_rows,col_filename}{:}];
                        if sum(row_content==mu_curently_considered_2) >= 1
                            similar_idx_2 = iterating_through_table_rows;
                            break;
                        else
                            continue;
                        end
                    end
                    % Detect if there is a confict
                    if similar_idx_1 == similar_idx_2
                        similar_idx = similar_idx_1;
                    else
                        warning(strcat("CONFLICT FOUND !",...
                            " MU#",num2str(mu_curently_considered_1),...
                            " of file ",row_filename,...
                            " and MU#",num2str(mu_curently_considered_2),...
                            " of file ",col_filename));
                        similar_idx = similar_idx_1;
                    end
                    MUs_matched{similar_idx,row_filename}= {[MUs_matched{similar_idx,row_filename}{:},mu_curently_considered_1]};
                    MUs_matched{similar_idx,col_filename} = {[MUs_matched{similar_idx,col_filename}{:},mu_curently_considered_2]};
                    MUs_matched{similar_idx,"Grid"}{:}(end+1:end+2) = [mu_grid_1,mu_grid_2];

                end
            end

        end
    end
end

%Keep only the list of unique grid indexes
MUs_matched{:,"Grid"} = cellfun(@unique,MUs_matched{:,"Grid"},'UniformOutput',false);
for mui=1:size(MUs_matched,1)
    if isempty(MUs_matched{mui,"Grid"}{:})
        MUs_matched{mui,"Muscle"} = {[]};
    elseif sum(MUs_matched{mui,"Grid"}{:}<=4) >= numel(MUs_matched{mui,"Grid"}{:})
        MUs_matched{mui,"Muscle"} = {'VL'};
    elseif sum(MUs_matched{mui,"Grid"}{:}>=5) >= numel(MUs_matched{mui,"Grid"}{:})
        MUs_matched{mui,"Muscle"} = {'VM'};
    else
        MUs_matched{mui,"Muscle"} = {'Conflict'};
    end
end
%Keep only the unique MUidx in each file column
for filei=1:size(files,1)
    MUs_matched{:,files(filei,2)} = cellfun(@unique,MUs_matched{:,files(filei,2)},'UniformOutput',false);
end

for mui=1:size(MUs_matched,1)
    mu_found_in_how_many_files = 0;
    for filei=1:size(files,1)
        if [MUs_matched{mui,files{filei,2}}{:}] > 0
            mu_found_in_how_many_files = mu_found_in_how_many_files + 1;
        end
    end
    MUs_matched{mui,"MU_found_in_how_many_files"} = {mu_found_in_how_many_files};
end

% Assign an arbitrary index to each matched MU
% Error when there is more than one element in a cell !
max_nb_of_matched_MUs = max_nb_of_MU;
for mui=1:size(MUs_matched,1)
    if isempty(MUs_matched{mui,"Muscle"}{:})
        max_nb_of_matched_MUs = mui-1;
        break;
    end
end
MUs_matched{1:max_nb_of_matched_MUs,'MU_idx'} = num2cell(1:max_nb_of_matched_MUs)';

%% SAVING OUTPUT
mkdir(output_dir_name)
save(strcat(output_dir_name,'\','MATCH_output.mat'), vars_to_save{:});

%% PLOT RESULTS
x = 1:size(files,1)+1; % +1 for total at the end
y_per_file = zeros(2,size(files,1)); %1st dim is number of MU per file (or total) ; 2nd dim is number of MU matched at least once
y_all_files = zeros(2,size(files,1)-1);
y_all_files(1,:) = 2:size(y_all_files,2)+1;
for total_yi = 1:size(y_all_files,2)
    y_all_files(2,total_yi) = sum([MUs_matched{:,"MU_found_in_how_many_files"}{:}]...
        >=y_all_files(1,total_yi));
end
for xi=1:numel(x)
    if xi<numel(x)
        %xlab{xi} = strrep(files{xi,2}(1:end-4),"_"," ");
        xlab{xi} = strrep(files{xi,2}(1:end),"_"," ");
        % need to find the matching column in "nb_of_MU_ach_file"
        y_per_file(1,xi) = nb_of_MU_each_file{2,xi}; % total nb of MUs in the file
        y_per_file(2,xi) = sum([MUs_matched{:,files{xi,2}}{:}]>0); % MUs matched at least once with the other files
    else
        xlab{xi} = "TOTAL";
    end
end

figure(100)
hold off
barcolors = zeros(size(files,1),3);
for bari=1:numel(x)
    if bari < numel(x)
        barcolors(bari,:) = [.5,.7,1];
    else
        barcolors(bari,:) = [1,.5,.3];
    end
end
barplot_handle = bar(x(1:end-1),y_per_file(1,:),'FaceColor',barcolors(1,:),...
    'FaceAlpha',0.5);
barplot_handle.CData = barcolors;
hold on
bar(x(1:end-1),y_per_file(2,:),'FaceColor',barcolors(1,:),...
    'FaceAlpha',0.8);
legend_labels = {"Total MUs in file";"MUs in file matched at least once"};
for total_yi = 1:size(y_all_files,2)
    bar(x(end),y_all_files(2,total_yi),'FaceColor',barcolors(end,:),...
        'FaceAlpha',(total_yi/size(y_all_files,2)) );
    if total_yi < size(y_all_files,2)
        legend_labels(end+1) = {strcat("MUs matched ", num2str(y_all_files(1,total_yi)) ,...
            " times or more")};
    else
         legend_labels(end+1) = {strcat("MUs matched", num2str(y_all_files(1,total_yi)) ,...
             " times (found in all files)")};
    end
end
legend(legend_labels);
ylabel("MU number");
xticks(1:numel(x));
xtickangle(5);
xticklabels(xlab);
xLim = xlim;
%line(xLim,[max(y_per_file(1,:)) max(y_per_file(1,:))],'Linew',2,'Linestyle','--','Color',[0,.3,.7]);
%line(xLim,[y_per_file(1,end) y_per_file(1,end)],'Linew',2,'Linestyle','--','Color',[.9,.1,0]);

% SAVE FIGURES
savefig(strcat(output_dir_name,'\','MATCH_output'));
saveas(figure(100),strcat(output_dir_name,'\','MATCH_output.png'));


% Generate new files (duplicate each file, keeping only the MU present in
% all files, and removing the other MUs)
% This needs to be done later, before that, creating figures to show what
% it looks like

%% CHECK MUAPS OF MATCHED MUs FOR EACH FILE

%close all
MUAP_shapes_to_save = cell(1,size(files,1));

if plot_MUAPs
    MUAP_window_size = round(0.05*2048);
    %MUAP_window_size = round(0.025*2048);

    condition_names = MUs_matched.Properties.VariableNames;
    condition_names(end) = [];
    condition_names(1:3) = [];

    for filei=1:size(files,1)
        eval( [condition_names{filei},' = load([files{filei,1},files{filei,2}]);'] );

        signal_file = eval([condition_names{filei},'.signal']);
        parameters_file = eval([condition_names{filei},'.parameters']);
        edition_file = eval([condition_names{filei},'.edition']);
        
        for gridi = 1:signal_file.ngrid
             SIGNAL = {}; % reset signal variable
            for channeli = 1:parameters_file.nbelectrodes
                SIGNAL{signal_file.coordinates{gridi}(channeli,1), signal_file.coordinates{gridi}(channeli,2)} = signal_file.data((gridi-1)*parameters_file.nbelectrodes+channeli,:);
                discardChannelsVec(signal_file.coordinates{gridi}(channeli,1), signal_file.coordinates{gridi}(channeli,2)) = signal_file.EMGmask{gridi}(channeli);
            end
            %SIGNAL = flip(SIGNAL);

            % Discard channels BEFORE doing differential, and for the
            % differential, check with next to next MU instead of just the
            % next one (if NaN)
            % Discard channels to be discarded (set to nan)
            for ch=1:numel(discardChannelsVec)
                if discardChannelsVec(ch)
                    %SIGNAL{ch} = [];
                    SIGNAL{ch} = NaN(1,numel(SIGNAL{ch}));
                end
            end

            SIGNAL_differential = cell(12,5);
            for coli=1:size(SIGNAL,2)
                temp_diff_col = cell2mat(SIGNAL(:,coli));
                temp_diff_col = diff(temp_diff_col);
                SIGNAL_differential(1:size(temp_diff_col,1),coli) = mat2cell(temp_diff_col, ones(size(temp_diff_col,1), 1), size(temp_diff_col,2));
            end
            clearvars coli temp_diff_col

            nb_of_MUs_in_grid = size(edition_file.Pulsetrain{gridi},1);
            for mui=1:nb_of_MUs_in_grid

                MUAP_shape{gridi,mui} = [];
                MUAPs_temp_mean = [];
                Peak_to_peak_temp = [];

                ch = 1;
                for electrode_grid_row=1:size(SIGNAL_differential,1)
                    for electrode_grid_col=1:size(SIGNAL_differential,2)
                        if ~isempty(SIGNAL_differential{electrode_grid_row,electrode_grid_col})
                            MUAP_shape_temp = check_MUAP(edition_file.Dischargetimes{gridi,mui},MUAP_window_size,SIGNAL_differential{electrode_grid_row,electrode_grid_col});
                            MUAPs_temp_mean{1}(ch,:) = mean(MUAP_shape_temp,1);
                            Peak_to_peak_temp{1}(ch) = max(MUAPs_temp_mean{1}(ch,:)) - min(MUAPs_temp_mean{1}(ch,:));
                            ch = ch+1;
                        end
                    end
                end
                idxmax = find(Peak_to_peak_temp{1}==max(Peak_to_peak_temp{1}(:)));
                [x,y] = max(MUAPs_temp_mean{1}(idxmax,:));
                ch = 1;
                MUAP_shape{gridi,mui} = cell(12,5);
                for electrode_grid_row=1:size(SIGNAL_differential,1)
                    for electrode_grid_col=1:size(SIGNAL_differential,2)
                        if ~isempty(SIGNAL_differential{electrode_grid_row,electrode_grid_col})
                            % Setting the MUAP shape time window that will
                            % be displayed
                            %MUAP_shape{gridi,mui}{electrode_grid_row,electrode_grid_col} = MUAPs_temp_mean{1}(ch,max([1,y-25]):min([y+25,size(MUAPs_temp_mean{1},2)]) );
                            %MUAP_shape{gridi,mui}{electrode_grid_row,electrode_grid_col} = MUAPs_temp_mean{1}(ch,max([1,y-50]):min([y+50,size(MUAPs_temp_mean{1},2)]) );
                            MUAP_shape{gridi,mui}{electrode_grid_row,electrode_grid_col} = MUAPs_temp_mean{1}(ch,1:size(MUAPs_temp_mean{1},2) );
                            ch = ch+1;
                        end
                    end
                end

            end

        end
        eval(strcat([condition_names{filei}],'.edition.MUAP_shapes = MUAP_shape;'));
        MUAP_shapes_to_save{filei} = MUAP_shape;
    end

% Save MUAP shapes
% re-calculate window size to make it the actual size used
    %MUAP_window_size = size(MUAP_shape{1}{1},2); % total window size, should be an odd number
MUAP_window_size = MUAP_window_size*2+1;
save(strcat(output_dir_name,'\','MUAP_shapes.mat'), plot_vars_to_save{:});

end

%% PLOT MUAPS SHAPES FOR EACH MU, AND SEE HOW DIFFERNET THEY ARE ACCORDING TO THE FILE

%close all

if plot_MUAPs
    MUAP_shape_dir_figures_output = "MUAP_shapes_figures";
    mkdir(strcat(output_dir_name,'\',MUAP_shape_dir_figures_output));

    MUAP_shapes_to_plot = {};
    colors_per_file = lines(size(files,1));
    for mui=1:max_nb_of_matched_MUs % reduce here for testing purposes
        if ~isempty(MUs_matched{mui,'MU_idx'}{:})
            MUAP_shapes_to_plot = cell(size(files,1),1);
            yMinMax = [-10,10]; % assign a default value in case there are empty values
            for filei=1:size(MUAP_shapes_to_plot,1)
                if ~isempty(MUs_matched{mui,condition_names{filei}}{:})
                    MUAP_shapes_to_plot{filei} = eval([condition_names{filei},'.edition.MUAP_shapes{'...
                        ,num2str(floor(MUs_matched{mui,condition_names{filei}}{:}/100)),',',num2str(mod(MUs_matched{mui,condition_names{filei}}{:},100)),'};']);
                    tempMin = cellfun(@min,[MUAP_shapes_to_plot{:}],'UniformOutput',false);
                    tempMax = cellfun(@max,[MUAP_shapes_to_plot{:}],'UniformOutput',false);
                    yMinMax = [ min(yMinMax(1),floor(min([tempMin{:}]))) , max(yMinMax(2),ceil(max([tempMax{:}]))) ];
                end
            end
            figure(mui)
            sgtitle(strcat("MU #",num2str(mui)," ; from ",MUs_matched{mui,"Muscle"}{:},...
                " ; grid(s) # [", num2str([MUs_matched{mui,"Grid"}{:}]),"]"));
        end
        legend_displayed = false;
        for filei=1:size(files,1)
            ch = 1;
            for rowi=1:size(MUAP_shapes_to_plot{filei},1)
                for coli=1:size(MUAP_shapes_to_plot{filei},2)
                    % Assign subplot
                    subplot(size(MUAP_shapes_to_plot{filei},1),size(MUAP_shapes_to_plot{filei},2),ch)
                    hold on
                    % Actually plot MUAP shape in subplot
                    plot(1:size(MUAP_shapes_to_plot{filei}{rowi,coli},2),MUAP_shapes_to_plot{filei}{rowi,coli},...
                        'LineWidth',(size(files,1)+1)-(filei/1.5),'color',colors_per_file(filei,:));
                    ch = ch+1;
                    ylim([yMinMax(1),yMinMax(2)]);
                    yticks(round(linspace(yMinMax(1), yMinMax(2), 4)));
                    ylabel('\mu V');
                    xticks([]);
                    
                    % Complicated thing to display legend. VERY INEFFICIENT
                    % Display legend only if already went through all files
                    non_empty_MUAP = find(~cellfun(@isempty,MUAP_shapes_to_plot));
                    last_non_empty_MUAP = non_empty_MUAP(end);
                    if filei==last_non_empty_MUAP && ~legend_displayed
                        legendLabels = {};
                        % Assign legend labels
                        for legend_entry_i = 1:numel(non_empty_MUAP)
                            legend_entry_corresponding_index = non_empty_MUAP(legend_entry_i);
                            legendLabels{legend_entry_i} = strcat(...
                                "MUAP shape of MU#",num2str(mod(MUs_matched{mui,condition_names{legend_entry_corresponding_index}}{:},100)),...
                                " in grid#",num2str(floor(MUs_matched{mui,condition_names{legend_entry_corresponding_index}}{:}/100)),...
                                " from file ",strrep(condition_names{legend_entry_corresponding_index},"_"," "));
                        end
                        % Display legend if legend is not already displayed
                        if ~legend_displayed
                            plotLegend = legend(legendLabels);
                        end
                        % If the legend displayed doesn't display all the
                        % MUAP shapes (not the same channels selected for
                        % the files), then delete the legend (and will try
                        % again in a different subplot)
                        % If all MUAP shapes displayed, set
                        % "legend_displayed" to true and stop trying a new
                        % subplot for which to display legend
                        if numel(plotLegend.String)==sum(~cellfun(@isempty, MUAP_shapes_to_plot))
                            legend_displayed = true;
                            plotLegend.Position(1:2) = [-0.2,0.9];
                        else
                            legend_displayed = false;
                            delete(plotLegend);
                        end
                    end

                end % end of "coli" loop
            end % end of "rowi" loop
        end % end of "filei" loop
                % FOR DEBUGGINF PURPOSE
                %             if mui==29
                %                 disp("STOP");
                %             end
    % SAVE FIGURES
    savefig(strcat(output_dir_name,'\',MUAP_shape_dir_figures_output,'\','MUAP_shapes_of_matched_MU_',num2str(mui)));
    saveas(figure(mui),strcat(output_dir_name,'\',MUAP_shape_dir_figures_output,'\','MUAP_shapes_of_matched_MU_',num2str(mui),'.png'));
    end % end of "muii" loop
end %% en of "if plot_MUAP"
