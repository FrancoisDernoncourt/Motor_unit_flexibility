clear all
clc
[filenames, path] = ...
    uigetfile(".mat",strcat("Select all online files to convert"), ...
    "Multiselect","on");

[filename_param,path_param] = ...
    uigetfile(".mat",strcat("Select edited online file with the 'parameters' struct"), ...
    "Multiselect","off");

loaded_file = strcat(path_param,filename_param);
load(loaded_file);
original_EMG_mask = signal.EMGmask;
vars_to_save = {"signal"}; %,"parameters"};

% Add parameters
if ~iscell(filenames)
    num_iter = 1;
else
    num_iter = numel(filenames);
end
for filei=1:num_iter
    if ~iscell(filenames)
        loaded_file = strcat(path,filenames);
    else
        loaded_file = strcat(path,filenames{filei});
    end
    disp(strcat("(",num2str(filei),"/",num2str(num_iter),") loading ",loaded_file));
    load(loaded_file);
%     signal.Dischargetimes = cellfun(@(x) x', signal.Dischargetimes, UniformOutput=false);
%     signal.emgtype = [1,1,1,1];
%     disp(strcat("(",num2str(filei),"/",num2str(num_iter),") saving ",loaded_file));
    signal.EMGmask = original_EMG_mask;
    save(loaded_file,vars_to_save{:});
end
