close all
clear all

file = uigetfile('Select edited file (Avrillon format ; duplicated removed)',"MultiSelect","off");
load(file);
%files_saved_suffix = "HANNING";
files_saved_suffix = "_windowed";


% Does the file contains the filteed signal? (toehr script of step 3_4 has
% already been used or not)
file_contains_filtered_signal = false;

window_already_selected = false;
if window_already_selected
    limOfEachWindow = [32762	74093
        157090	198085
        282090	320397
        402050	443380];
end

% force_file = uigetfile('Select the signal containing the force file',"MultiSelect","off");
% load(force_file);

% removeDuplicatedOutput_file = uigetfile('Select the file with removed duplicate output',"MultiSelect","off");
% load(removeDuplicatesOutput_file);

%% PARAMETERS FOR SELECTION
samplingRate = 2048;
clickSelectOrWindowSelect = true; %true for click select = choosing manually begining and end of the window
windowSize = 10; %in s ; used only if "clickSelectOrWindowSelect = false"
%figureName = "S1_hufr_plateau5_IDR_LowPass";
figureName = "window_cut";
alreadyBinarized = false;
Predetermined_cuts = []; %if ~= [] ; then cuts will already be made. Works only if clickSelectOrWindowSelect = true
signalDuration = length(signal.path);
% originalSavename = "P4_CONT"; %savename(1:end-16);
originalSavename = file(1:end-4);

%% DISPLAY SPIKE TRAINS AND FORCE SIGNAL

if ~window_already_selected

    edition.binaryDischargeMatrix = cell(1,signal.ngrid);
    mu_count = 0;
    for i=1:signal.ngrid
        edition.binaryDischargeMatrix{i} = zeros(size(edition.Pulsetrain{i}));
        for mui=1:size(edition.Pulsetrain{i},1)
            for spike_idx=1:size(edition.Dischargetimes{i,mui},2)
                spike_time = edition.Dischargetimes{i,mui}(spike_idx);
                edition.binaryDischargeMatrix{i}(mui,spike_time) = 1;
            end
            mu_count = mu_count+1;
            concat_binary_discharge_matrix(mu_count,:) = edition.binaryDischargeMatrix{i}(mui,:);
            which_MU_belongs_to_which_grid(mu_count) = i;
        end
    end

    %plot_handle = plot(signal.path.*length(which_MU_belongs_to_which_grid),"LineWidth",3,'Color',[.5,.5,.5]);
    hold on
    grid_muscles_correspondance = [1,nan,nan,nan;2,nan,nan,nan;3,4,5,6];
    muscles_colors = lines(size(grid_muscles_correspondance,1));
    mu_count = 0;
    for i=1:signal.ngrid
        [which_muscle,~] = find(i==grid_muscles_correspondance);
        color_to_use = muscles_colors(which_muscle,:);
        % color_to_use(4) = 0.5; % I tried to have alpha data but it didn't work
        for mui=1:size(edition.binaryDischargeMatrix{i},1)
            vec_to_plot = edition.binaryDischargeMatrix{i}(mui,:);
            vec_to_plot(vec_to_plot==0) = NaN;
            vec_to_plot = vec_to_plot.* mu_count;

            mu_count = mu_count + 1;
            plot(1:length(signal.path),vec_to_plot,'|','Color',color_to_use);
        end
    end
%     % force instead of path = 64 * number of grids + 1
%     plot((signal.path ./ max(signal.path))...
%         .*length(which_MU_belongs_to_which_grid),"LineWidth",1.5,'Color',[0,0,0]);
      plot((signal.data(64*signal.ngrid+1,:) ./ min(signal.data(64*signal.ngrid+1,:)) )...
        .*length(which_MU_belongs_to_which_grid),"LineWidth",1.5,'Color',[0,0,0]);

    title({['Please select the windows IN CHONOLOGICAL ORDER. Press [SPACE] when selection is finished']})
    xTicksTemp = xticks;
    xTicksLabelsTemp = xticks ./ 2048; %Change the x axis to show seconds instead of time points
    xticklabels((round(xTicksLabelsTemp*10))/10) %round the x values to 1 digit after the point
    xlabel(['Time in seconds'])
    ylabel('Motor units')

    xLim = xlim;
    yLim = ylim;


    %% SELECT WINDOWS

    selecting = true;
    inputButton = nan;
    createPatch = false;
    currentCut = nan;
    cutXcoord = nan;
    if exist("patches","var")
        clearvars patches
    end

    set(gcf,"KeyPressFcn", @KeyDownCallback);
    set(gcf,"WindowButtonDownFcn", @(hObject, eventdata, selectOption, winSize, fsamp) ...
        MouseClickCallback(hObject, eventdata, clickSelectOrWindowSelect, windowSize, samplingRate) );

    if ~clickSelectOrWindowSelect
        cursorPatch = patch([0 (windowSize*samplingRate) (windowSize*samplingRate) 0],[yLim(1) yLim(1) yLim(2) yLim(2)],'c'); %Draw the initial patch outside the window
        cursorPatch.FaceVertexAlphaData = 0.1;
        cursorPatch.FaceAlpha = 'flat';
        cursorPatch.EdgeAlpha = 0;
        % Update patch position according to cursor
        set(gcf, 'WindowButtonMotionFcn', @(object, eventdata, patch, size) mouseMove(object, eventdata, cursorPatch, windowSize*samplingRate) );
    end

    if ~isempty(Predetermined_cuts)
        cutsExecuted=1;
    end

    while (selecting)

        if ~isempty(Predetermined_cuts) && cutsExecuted <= length(Predetermined_cuts)
            assignedCutsCallback(Predetermined_cuts(cutsExecuted))
            cutsExecuted=cutsExecuted+1;
        end

        % Space is pressed and selection is finished
        if inputButton == 'space'
            selecting = false;
        end

        % Create patches
        if createPatch

            if isnan(cutXcoord)
                cutXcoord = [];
            end

            for x=1:length(currentCut)
                cutXcoord(end+1) = currentCut(x);
            end

            % Check if the cut was made in chronological order
            if length(cutXcoord) > 1
                if cutXcoord(end) < cutXcoord(end-1)
                    error("Please select the windows in CHRONOLOGICAL order");
                end
            end

            if clickSelectOrWindowSelect
                if length(cutXcoord) == 1
                    patches(1) = patch([0 xLim(1) xLim(2) 0],[yLim(1) yLim(1) yLim(2) yLim(2)],'k'); % Draw the initial patches outside the window
                elseif mod(length(cutXcoord),2) ~= 1 % if nb of cuts > 1, add a patch only when even number of cuts
                    patches(ceil((length(cutXcoord)/2)+0.1)) = patch([0 xLim(1) xLim(2) 0],[yLim(1) yLim(1) yLim(2) yLim(2)],'k'); % Draw the initial patches outside the window
                end
            else
                if length(cutXcoord) == 2
                    patches(1) = patch([0 xLim(1) xLim(2) 0],[yLim(1) yLim(1) yLim(2) yLim(2)],'k');
                    patches(2) = patch([0 xLim(1) xLim(2) 0],[yLim(1) yLim(1) yLim(2) yLim(2)],'k');
                else %if one window already exists, one patch is added per click
                    patches(length(patches)+1) = patch([0 xLim(1) xLim(2) 0],[yLim(1) yLim(1) yLim(2) yLim(2)],'k');
                end
            end
            createPatch = false;

            % Draw patches
            for i=1:length(patches)
                % Each patch xlim(1)
                if (2*(i-1)) < 1 % 1st patch has x value starting at 0
                    patchMinLim = 1; % 1 and not zero so that the array is not out of range when deleting
                else
                    patchMinLim = cutXcoord((2*(i-1))); %(2*(i-1))=2 mais cutXcoord n'a qu'une seule valeur pour l'instant !
                end
                % Each patch xlim(2)
                if (2*(i-1)) >= length(cutXcoord) % last patch has x value ending at xlim
                    patchMaxLim = xLim(2);
                else
                    patchMaxLim = cutXcoord((2*(i-1))+1);
                end

                set(patches(i),'XData',[patchMinLim, patchMaxLim , patchMaxLim , patchMinLim]);
                patches(i).FaceVertexAlphaData = 0.1;
                patches(i).FaceAlpha = 'flat';
                patches(i).EdgeAlpha = 0;
            end

        end
        pause(0.2);
    end

    clearvars MUPulsesTotal

    %% SAVE INDIVIDUAL WINDOWS
    if ~isempty(cutXcoord)
        numberOfWindows = ceil(length(cutXcoord)/2);
    else
        numberOfWindows = 1;
    end

    limOfEachWindow = zeros(numberOfWindows,2);

    if isempty(patches)
        limOfEachWindow(1,1) = 1;
        limOfEachWindow(1,2) = signalDuration;
    else
        for p=1:length(patches) % "p" for "patch"

            if length(patches) == numberOfWindows
                limOfEachWindow(p,1) = patches(p).XData(2);
                if p == length(patches)
                    limOfEachWindow(p,2) = signalDuration;
                else
                    limOfEachWindow(p,2) = patches(p+1).XData(1);
                end
            elseif length(patches) > numberOfWindows
                if p <= (length(patches)-1)
                    limOfEachWindow(p,1) = patches(p).XData(2);
                    limOfEachWindow(p,2) = patches(p+1).XData(1);
                else
                    break
                end
            end

        end
    end

    save(strcat(originalSavename,"_selected_windows_timepoints",files_saved_suffix,".mat"),"limOfEachWindow");

end

% %% CUT INDIVIDUAL WINDOWS AND SAVE FILES
% 
% MUindex = 1;
% binarizedMUPulses = cell(1, signal.ngrid);
% if file_contains_filtered_signal
%     filteredMUs = cell(1, signal.ngrid );
% end
% MUsIPTs = cell(1, signal.ngrid );
% for i=1:signal.ngrid
%     binarizedMUPulses{i} = zeros(size(edition.Pulsetrain{i},1), signalDuration);
%     filteredMUs{i} = zeros(size(edition.Pulsetrain{i},1), signalDuration );
%     MUsIPTs{i} = zeros(size(edition.Pulsetrain{i},1), signalDuration );
%     for j=1:size(binarizedMUPulses{i},1)
%         % Binarized spike vectors of each MUs
%         binarizedMUPulses{i}(j,:) = edition.binaryDischargeMatrix{i}(j,:);
%         % filtered signal of each MU
%         if file_contains_filtered_signal
%             filteredMUs{i}(j,:) = edition.filteredMUs{i}(j,:);
%         end
%         % IPTs, Pulsetrain
%         MUsIPTs{i}(j,:) = edition.Pulsetrain{i}(j,:);
%         MUindex = MUindex+1; % actually not used
%     end
% end
% % Dischargetimes of each MUs
% dischargeTimesToCut = edition.Dischargetimes;
% % add something for signal.data
% emgSignalData = signal.data;
% % Path and target
% pathToCut = signal.path;
% targetToCut = signal.target;
% force = signal.path;
% 
% if window_already_selected
%     numberOfWindows = size(limOfEachWindow,1);
% end
% 
% selectedWindows = cell(1,numberOfWindows);
% for w=1:numberOfWindows % "w" for "window"
%     for g=1:signal.ngrid % "g" for "grid"
%         % Binary spike trains
%         temp = [];
%         for MUi = 1:size(MUsIPTs{g},1)
%             temp(MUi,:) = binarizedMUPulses{g}(MUi,limOfEachWindow(w,1):limOfEachWindow(w,2));
%         end
%         selectedWindows{w}.edition.binarySpikeTrains{g} = temp;
% 
%         % Recalculated spike discharge times
%         %         temp = dischargeTimesToCut;
%         % %         temp = cell(size(dischargeTimesToCut));
%         % %         for MUi = 1:size(MUsIPTs{g},1) %Maybe could be simplified by using cellfun
%         % %             temp{g,MUi} = dischargeTimesToCut{g,MUi};
%         % %             temp{g,MUi} = temp{g,MUi}(temp{MUi}>limOfEachWindow(w,1));
%         % %             temp{g,MUi} = temp{g,MUi}(temp{MUi}<limOfEachWindow(w,2));
%         % %             temp{g,MUi} = temp{g,MUi} - limOfEachWindow(w,1);
%         % %         end
%         %         selectedWindows{w}.edition.Dischargetimes{g,:} = temp{g,:};
% 
%         % Filtered MU signals
%         if file_contains_filtered_signal
%             temp = [];
%             for MUi = 1:size(MUsIPTs{g},1)
%                 temp(MUi,:) = filteredMUs{g}(MUi,limOfEachWindow(w,1):limOfEachWindow(w,2));
%             end
%             selectedWindows{w}.edition.filteredMUs{g} = temp;
%         end
% 
%         % IPTs
%         temp = [];
%         for MUi = 1:size(MUsIPTs{g},1)
%             temp(MUi,:) = MUsIPTs{g}(MUi,limOfEachWindow(w,1):limOfEachWindow(w,2));
%         end
%         selectedWindows{w}.edition.Pulsetrain{g} = temp;
%     end
%     % Discharge times
%     temp = [];
%     temp = dischargeTimesToCut;
%     cellfx = @(x) x(x>limOfEachWindow(w,1) & x<limOfEachWindow(w,2));
%     temp = cellfun(cellfx, temp, 'UniformOutput', false);
%     cellfx = @(x) x-limOfEachWindow(w,1);
%     temp = cellfun(cellfx, temp, 'UniformOutput', false);
%     selectedWindows{w}.edition.Dischargetimes = temp;
%     clearvars cellfx
%     % EMG data
%     selectedWindows{w}.signal.data = emgSignalData(:,limOfEachWindow(w,1):limOfEachWindow(w,2));
%     % Path and target
%     selectedWindows{w}.signal.target = pathToCut(limOfEachWindow(w,1):limOfEachWindow(w,2));
%     selectedWindows{w}.signal.path = targetToCut(limOfEachWindow(w,1):limOfEachWindow(w,2));
%     selectedWindows{w}.signal.path = force(limOfEachWindow(w,1):limOfEachWindow(w,2));
%     % Other values
%     selectedWindows{w}.signal.ngrid = signal.ngrid;
%     selectedWindows{w}.signal.nChan = signal.nChan;
%     selectedWindows{w}.signal.EMGmask = signal.EMGmask;
%     selectedWindows{w}.signal.emgtype = signal.emgtype;
%     selectedWindows{w}.signal.IED = signal.IED;
%     selectedWindows{w}.signal.coordinates = signal.coordinates;
%     selectedWindows{w}.signal.gridname = signal.gridname;
%     selectedWindows{w}.signal.muscle = signal.muscle;
%     selectedWindows{w}.signal.fsamp = signal.fsamp;
%     selectedWindows{w}.savename = strcat(originalSavename,'_window_',num2str(w));
%     % Save each window in a separate file
%     signal = selectedWindows{w}.signal;
%     edition = selectedWindows{w}.edition;
%     savename = selectedWindows{w}.savename;
%     % "parameters" is strictly the same for every file
%     save(strcat(selectedWindows{w}.savename,files_saved_suffix,".mat"), ...
%         "edition","savename","signal","parameters",'-v7.3'); %'-v7.3' makes sure everything gets saved even if the file is very big);
% end

%% FUNCTIONS CALLED

function mouseMove (~, ~, patch, size)
C = get (gca, 'CurrentPoint');
set(patch,'XData',[C(1,1), C(1,1) + size , C(1,1) + size , C(1,1)]);
end

function KeyDownCallback (~, keyPressedEventData)
coordinates = get(gca,'CurrentPoint');
xClickCoord = coordinates(1,1);
disp(xClickCoord);
assignin('base','xClickCoord',xClickCoord);
inputButton = keyPressedEventData.Key;
assignin('base','inputButton',inputButton);
end


function MouseClickCallback (~, ~, clickSelectOrWindowSelect, windowSize, samplingRate)

coordinates = get(gca,'CurrentPoint');
xClickCoord = coordinates(1,1);
%assignin('base','xClickCoord',xClickCoord);

cutXcoord(1) = round(xClickCoord);

% Make the equivalent of two cuts if window selection is on
if ~clickSelectOrWindowSelect
    cutXcoord(2) = round(xClickCoord)+(windowSize*samplingRate);
end

assignin('base','currentCut',cutXcoord);

createPatch = true;
assignin('base','createPatch',createPatch)

end

function assignedCutsCallback (coord)
assignin('base','currentCut',coord);
createPatch = true;
assignin('base','createPatch',createPatch)
end
