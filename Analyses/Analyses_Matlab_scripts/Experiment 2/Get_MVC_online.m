clc

[filename,path] = ...
    uigetfile(".mat",strcat("Select MVC file"), ...
    "Multiselect","off");

loaded_file = strcat(path,filename);
load(loaded_file);
vars_to_save = {"force_MVC","force_baseline"};

force_reversed = false;
if force_reversed
    Force = Force.*(-1);
end

figure()
plot(Force)
title("Please select the window edges over which the MVC force value will be determined");
[input_x, ~] = ginput(2);
input_x = sort(input_x,'ascend');
input_x = round(input_x);

% if force_reversed
%     force_MVC = min(Force(input_x(1):input_x(2))) - Offset_force;
% else
%     force_MVC = max(Force(input_x(1):input_x(2))) + Offset_force;
% end
force_MVC = max(Force(input_x(1):input_x(2))) + Offset_force;

title("Please select the window edges over which the baseline force will be computed");
[input_x, ~] = ginput(2);
input_x = sort(input_x,'ascend');
input_x = round(input_x);

% if force_reversed
%     force_baseline = mean(Force(input_x(1):input_x(2))) - Offset_force;
% else
%     force_baseline = mean(Force(input_x(1):input_x(2))) + Offset_force;
% end
force_baseline = mean(Force(input_x(1):input_x(2))) + Offset_force;

close(gcf());

save("Force_values",vars_to_save{:});