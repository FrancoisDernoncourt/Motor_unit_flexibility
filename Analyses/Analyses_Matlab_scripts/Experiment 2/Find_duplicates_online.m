% PSEUDO-CODE :
%   Select files
%   Assign a value to each MU according to which grid they belong to
%   Concatenate everything
%   remduplicates
%   Re-separate the files
%   save the generated files in a new folder

%% SELECT FILES
[filename,~] = uigetfile('*.mat', 'Select edited file (Avrillon format)', 'MultiSelect', 'off');
load(filename);
nbOfMUs = 0;
signal.ngrid = size(signal.Pulsetrain,2);
for i=1:signal.ngrid
    nbOfMUs = nbOfMUs + size(edition.Pulsetrain{1},1);
end
signalLength = length(signal.target);

%% CONCATENATE EVERYTHING
concat.MUPulses = cell(1,nbOfMUs);
concat.IPTs = zeros(nbOfMUs,signalLength);
concat.whichMUbelongsToWhichGrid = zeros(nbOfMUs,1);
currentMUiteratedOn = 0;
corresponding_MUs = [];
for i=1:signal.ngrid
    nbOfMUsInGrid = size(edition.Pulsetrain{i},1);
    for n=1:nbOfMUsInGrid
        currentMUiteratedOn = currentMUiteratedOn + 1;
        concat.MUPulses{currentMUiteratedOn} = edition.Dischargetimes{i,n};
        concat.IPTs(currentMUiteratedOn,:) = edition.Pulsetrain{i}(n,:);
        concat.whichMUbelongsToWhichGrid(currentMUiteratedOn) = i;
        corresponding_MUs(end+1) = i*100+n;
    end
end

%% REMOVE DUPILCATES
[IPTs_keep , MUPulses_keep, duplicatesList, removedList] = ...
    remove_duplicates(concat.IPTs, concat.MUPulses, concat.MUPulses, round(2048/10), 0.00025, 0.15, 2048);


%% FIND DUPLICATES
duplicatesList = duplicatesList';
highest_nb_of_duplicates = max(cellfun(@numel, duplicatesList));
if isempty(highest_nb_of_duplicates) || highest_nb_of_duplicates == 0
    highest_nb_of_duplicates = 1;
end
duplicates_table = nan(numel(duplicatesList),highest_nb_of_duplicates+2);
duplicate_iter = 0;
for mui=1:size(duplicates_table,1)
    duplicates_table(mui,1) = corresponding_MUs(mui);
    if isempty(duplicatesList{mui})
        continue
    end
    duplicates_table(mui, 2: 2+(numel(duplicatesList{mui})-1) ) = corresponding_MUs(...
        duplicatesList{mui});
    duplicate_iter = duplicate_iter + 1;
    if removedList{mui}==mui
        duplicates_table(mui,end) = duplicates_table(mui,2);
    else
        duplicates_table(mui,end) = duplicates_table(mui,1);
    end
end

if size(duplicates_table,2) > 3
    duplicates_table_temp = duplicates_table;
    for mui=1:size(duplicates_table,1)
        duplicates_table{mui,2} = [duplicates_table_temp{mui,2:end-1}];
    end
    duplicates_table(:,3:end-1) = [];
end

varNames = {'MU_idx','Duplicated_MU_idx','MU_idx_to_keep'};
duplicates_table = num2cell(duplicates_table);
duplicates_table = cell2table(duplicates_table,'VariableNames',varNames);

varsToSave = {"duplicates_table"};
save(strcat('Duplicates_list.mat'),varsToSave{:});
