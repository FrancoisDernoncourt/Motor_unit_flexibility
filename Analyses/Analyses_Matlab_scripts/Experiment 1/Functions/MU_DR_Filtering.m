function filteredBinaryDischargeMatrix = MU_DR_Filtering(BinaryDischargeMatrix,SamplingRate,isBinary,signalDuration,FilteringMethod,Edges,Normalization)
%This function filters the binary discharge matrix used as
%input and outputs a smoothed matrix correspond to smoothed spike trains of
%each MU
%   1st argument = matrix to filter (should be composed of zeros and 1)
%   2nd argument = sampling rate (in samples per s)
%   3rd argument = "true" if already a matrix in binary format ; "false" if cell array in non-binary format (just spike times)
%   4th argument (optional, necessary only if not already in binary format) = int value equal to the signal duration in samples
%   5th argument = filter method "IDR_LowPass" or "Hann_conv_HighPass" or "Hann_conv_Lowpass_Highpass" or "Gauss_conv". Default = "Hann_conv_HighPass"
%   6th argument = what to do with edges ? Options : "Remove" ; "Zero-pad"   pads with zeros ; "Nan-pad" pads with nans (Nan-pad is default option)
%   7th argument = normalization method. None by default. Options : "None ;   "MaxVal" ; "Between0and1" ; "Between-1and1" ; "max_IDR"

if ~isBinary
    if ~exist("signalDuration","var")
        error("if your input is non-binary, you must include your signal duration as an inpjt. Signal duration is an int value")
    else
        MUPulses = BinaryDischargeMatrix;
        BinaryDischargeMatrix = zeros(length(MUPulses),signalDuration);
        for i=1:size(BinaryDischargeMatrix,1)
            for spike=1:length(MUPulses{i})
                BinaryDischargeMatrix(i,MUPulses{i}(spike)) = 1;
            end
        end
    end
end

nbOfMUs = size(BinaryDischargeMatrix,1);
signalDuration = size(BinaryDischargeMatrix,2);

if FilteringMethod == "Hanning"

    % From DelVecchio code neural modules
    fsamp = SamplingRate; %set your fsamp;
    Wind_s = 1;  % hanning window duration 
    HanningW = 2/round(fsamp*Wind_s)*hann(round(fsamp*Wind_s)); %unitary area 

    filteredBinaryDischargeMatrix = nan(size(BinaryDischargeMatrix));
    for MUidx = 1:nbOfMUs
        Firings_temp = BinaryDischargeMatrix(MUidx,:);
        filteredBinaryDischargeMatrix(MUidx,:) = filtfilt(HanningW,1,Firings_temp*fsamp);
    end

elseif FilteringMethod == "IDR_LowPass"
    
    ISI_thresh_for_0_IDR = 1; %in seconds = amount of quiescent time for the MN above which the IDR is considered to be zero (derecruited MN)
    lowPassCutOff = 5; %in hz

    filteredBinaryDischargeMatrix = nan(size(BinaryDischargeMatrix));
    [b,a] = butter(1 , lowPassCutOff / (SamplingRate/2) , 'low' ); %1st order low-pass butterworth filter with cut-off frequency of 5Hz
    t = 1:signalDuration;

    % get instantaneous discharge rate
    for MUidx = 1:nbOfMUs
        spikesOfCurrentMU = find(BinaryDischargeMatrix(MUidx,:));
        IDRofSpikes = zeros(size(spikesOfCurrentMU));
        IDRofSpikes(2:end) = (1./diff(spikesOfCurrentMU))*SamplingRate;
        for IDRidx = 1:length(spikesOfCurrentMU)
            filteredBinaryDischargeMatrix(MUidx,spikesOfCurrentMU(IDRidx)) = IDRofSpikes(IDRidx);
        end
        % Record NaN values
        tempNanIndexes = find(isnan(interp1( t(~isnan(filteredBinaryDischargeMatrix(MUidx,:))) , filteredBinaryDischargeMatrix(MUidx,t(~isnan(filteredBinaryDischargeMatrix(MUidx,:)))) , t , 'linear')));
        last_valid_recruitment_idx = 1;
        first_valid_derecruitment_idx = signalDuration;
        % Interpolation
        % Interpolate only for chunks for which length(NaN) < ISI_thresh_for_0_IDR
%         distime = find(~isnan(filteredBinaryDischargeMatrix(MUidx,:)));
%         ISI(1) = nan;
%         ISI(2:length(distime)) = diff(distime);
        ISI = diff(spikesOfCurrentMU);
        last_spike_before_derecruitment_indexes = find(ISI > (ISI_thresh_for_0_IDR*SamplingRate));
        first_spike_after_derecruitment_indexes = last_spike_before_derecruitment_indexes+1;
        last_spike_before_derecruitment_times = spikesOfCurrentMU(last_spike_before_derecruitment_indexes);
        first_spike_after_derecruitment_times = spikesOfCurrentMU(first_spike_after_derecruitment_indexes);
        filteredBinaryDischargeMatrix(MUidx,:) = interp1( t(~isnan(filteredBinaryDischargeMatrix(MUidx,:))) , filteredBinaryDischargeMatrix(MUidx,t(~isnan(filteredBinaryDischargeMatrix(MUidx,:)))) , t , "linear" ,false); % false brings extrapolated values to 0 (NaN values)
        for derecruitmenti = 1:numel(last_spike_before_derecruitment_indexes)
            nb_of_derecruitment_samples = length(last_spike_before_derecruitment_times(derecruitmenti):first_spike_after_derecruitment_times(derecruitmenti));
            filteredBinaryDischargeMatrix(MUidx,last_spike_before_derecruitment_times(derecruitmenti):first_spike_after_derecruitment_times(derecruitmenti)) = zeros(1,nb_of_derecruitment_samples);
        end
        % Low-pass filter
        filteredBinaryDischargeMatrix(MUidx,:) = filtfilt(b,a,filteredBinaryDischargeMatrix(MUidx,:));
%         % High-pass filtered
%         filteredBinaryDischargeMatrix(MUidx,:) = filtfilt(d,c,filteredBinaryDischargeMatrix(MUidx,:));
        % Re-apply NaN values
        for nanIndexes=1:length(tempNanIndexes)
            if Edges == "" || Edges == "Nan_pad"
                filteredBinaryDischargeMatrix(MUidx,tempNanIndexes(nanIndexes)) = nan;
            elseif Edges == "Zero-pad"
                filteredBinaryDischargeMatrix(MUidx,tempNanIndexes(nanIndexes)) = 0;
            elseif Edges == "Remove"
                % find last valid recruitment index & first valid derecruitment index
                invalid_Idx = tempNanIndexes;
                idx_at_which_filtered_signal_is_valid = diff(invalid_Idx);
                idx_at_which_filtered_signal_is_valid = find(idx_at_which_filtered_signal_is_valid>1);
                invalid_Idx_before_Recruitment = invalid_Idx(1:idx_at_which_filtered_signal_is_valid);
                invalid_Idx_after_Derecruitment = invalid_Idx(invalid_Idx>idx_at_which_filtered_signal_is_valid);
                if first_valid_derecruitment_idx > min(invalid_Idx_after_Derecruitment)
                    first_valid_derecruitment_idx = min(invalid_Idx_after_Derecruitment);
                end
                if last_valid_recruitment_idx < max(invalid_Idx_before_Recruitment)
                    last_valid_recruitment_idx = max(invalid_Idx_before_Recruitment);
                end
            else
                error("Please assign a valid method to deal with edges")
            end
        end
    end

    if Edges == "Remove"
        if first_valid_derecruitment_idx < signalDuration
            filteredBinaryDischargeMatrix(:,end:-1:first_valid_derecruitment_idx) = [];
        end
        if last_valid_recruitment_idx > 1
            filteredBinaryDischargeMatrix(:,1:last_valid_recruitment_idx) = [];
        end
    end

elseif FilteringMethod == "Hann_conv_HighPass" || FilteringMethod == ""

    kernelDuration = 0.4; %in s
    highPassCutOff = 0.75; %in hz

    filteredBinaryDischargeMatrix = nan(size(BinaryDischargeMatrix)); %Needs the size to be changed
    samples_lost_to_conv = round((SamplingRate * kernelDuration)-1);
    kernel = hann(round(kernelDuration*SamplingRate));

    [b,a] = butter(3 , highPassCutOff / (SamplingRate/2) , 'high' ); % 3rd order high-pass butterworth filter with cut-off frequency of 0.75Hz

    for MUidx=1:nbOfMUs
        tempFilteredBinaryDischargeMatrix = conv(BinaryDischargeMatrix(MUidx,:),kernel,"valid");
        tempFilteredBinaryDischargeMatrix = filtfilt(b,a,tempFilteredBinaryDischargeMatrix);
        filteredBinaryDischargeMatrix(MUidx,round(samples_lost_to_conv/2)+1 : round(samples_lost_to_conv/2)+length(tempFilteredBinaryDischargeMatrix)) = tempFilteredBinaryDischargeMatrix;
    end

    % Padding options
    nan_matrix = isnan(filteredBinaryDischargeMatrix);
    if Edges == "" || Edges == "Nan_pad"
        % Do nothing ; already with NaN values
    elseif Edges == "Zero-pad"
        filteredBinaryDischargeMatrix(nan_matrix) = 0;
    elseif Edges == "Remove"
        nan_cols = nan_matrix(1,:);
        col_idx_to_remove = find(nan_cols);
        filteredBinaryDischargeMatrix(:,col_idx_to_remove) = [];
    else
        error("Please assign a valid method to deal with edges")
    end

elseif FilteringMethod == "Hann_conv_Lowpass_Highpass"

    kernelDuration = 0.4; %in s
    highPassCutOff = 0.75; %in hz
    lowPassCutOff = 5; %in hz

    filteredBinaryDischargeMatrix = nan(size(BinaryDischargeMatrix)); %Needs the size to be changed
    samples_lost_to_conv = round((SamplingRate * kernelDuration)-1);
    kernel = hann(round(kernelDuration*SamplingRate));

    [b,a] = butter(3 , highPassCutOff / (SamplingRate/2) , 'high' ); % 3rd order high-pass butterworth filter with cut-off frequency of 0.75Hz
    [d,c] = butter(1 , lowPassCutOff / (SamplingRate/2) , 'low' ); %1st order low-pass butterworth filter with cut-off frequency of 5Hz

    for MUidx=1:nbOfMUs
        tempFilteredBinaryDischargeMatrix = conv(BinaryDischargeMatrix(MUidx,:),kernel,"valid");
        tempFilteredBinaryDischargeMatrix = filtfilt(b,a,tempFilteredBinaryDischargeMatrix);
        tempFilteredBinaryDischargeMatrix = filtfilt(d,c,tempFilteredBinaryDischargeMatrix);
        filteredBinaryDischargeMatrix(MUidx,round(samples_lost_to_conv/2)+1 : round(samples_lost_to_conv/2)+length(tempFilteredBinaryDischargeMatrix)) = tempFilteredBinaryDischargeMatrix;
    end

    % Padding options
    nan_matrix = isnan(filteredBinaryDischargeMatrix);
    if Edges == "" || Edges == "Nan_pad"
        % Do nothing ; already with NaN values
    elseif Edges == "Zero-pad"
        filteredBinaryDischargeMatrix(nan_matrix) = 0;
    elseif Edges == "Remove"
        nan_cols = nan_matrix(1,:);
        col_idx_to_remove = find(nan_cols);
        filteredBinaryDischargeMatrix(:,col_idx_to_remove) = [];
    else
        error("Please assign a valid method to deal with edges")
    end

elseif FilteringMethod == "Gauss_conv"

    kernelDuration = 0.025; %in s

    filteredBinaryDischargeMatrix = nan(size(BinaryDischargeMatrix)); %Needs the size to be changed
    samples_lost_to_conv = round((SamplingRate * kernelDuration)-1);
    kernelSize = SamplingRate * kernelDuration;  % Size of the kernel (odd number)
    % Create a 1D array for the kernel indices
    indices = -(kernelSize-1)/2 : (kernelSize-1)/2;
    % Create the 1D Gaussian kernel using normpdf
    sigma = 1;
    kernel = normpdf(indices, 0, sigma);
    % Normalize the kernel to have a sum of 1
    kernel = kernel / sum(kernel);

    kernel = fspecial('gaussian', round(kernelDuration*SamplingRate));
    for MUidx=1:nbOfMUs
        tempFilteredBinaryDischargeMatrix = conv(BinaryDischargeMatrix(MUidx,:),kernel,"valid");
        filteredBinaryDischargeMatrix(MUidx,round(samples_lost_to_conv/2)+1 : round(samples_lost_to_conv/2)+length(tempFilteredBinaryDischargeMatrix)) = tempFilteredBinaryDischargeMatrix;
    end

    % Padding options
    nan_matrix = isnan(filteredBinaryDischargeMatrix);
    if Edges == "" || Edges == "Nan_pad"
        % Do nothing ; already with NaN values
    elseif Edges == "Zero-pad"
        filteredBinaryDischargeMatrix(nan_matrix) = 0;
    elseif Edges == "Remove"
        nan_cols = nan_matrix(1,:);
        col_idx_to_remove = find(nan_cols);
        filteredBinaryDischargeMatrix(:,col_idx_to_remove) = [];
    else
        error("Please assign a valid method to deal with edges")
    end

else
    error("Please select a valid filtering method")
end

if Normalization == "" || Normalization == "None"
else
    if Normalization == "MaxVal"
        for i=1:size(filteredBinaryDischargeMatrix,1)
            filteredBinaryDischargeMatrix(i,:) = filteredBinaryDischargeMatrix(i,:) ./ max(filteredBinaryDischargeMatrix(i,:));
        end
    elseif Normalization == "Between0and1"
        for i=1:size(filteredBinaryDischargeMatrix,1)
            warning("Chosen normalization method is not integrated yet");
        end
    elseif Normalization == "Between-1and1"
        for i=1:size(filteredBinaryDischargeMatrix,1)
            warning("Chosen normalization method is not integrated yet");
        end
    elseif Normalization == "Max_IDR"
        for i = 1:size(filteredBinaryDischargeMatrix,1)
            % determine max IDR
            nbOfDischargesTemp = unique(cumsum(BinaryDischargeMatrix(i,:)));
            IDRsTemp = histcounts(cumsum(BinaryDischargeMatrix(i,:)) , length(nbOfDischargesTemp));
            IDRsTemp = IDRsTemp(1:(end-1));
            IDRsTemp(1) = NaN;
            max_IDR = 2048/min(IDRsTemp);
            min_IDR = 2048/max(IDRsTemp);
            % Normalize filtered signal according to max IDR
            filteredBinaryDischargeMatrix(i,:) = rescale(filteredBinaryDischargeMatrix(i,:),min_IDR,max_IDR);
            %filteredBinaryDischargeMatrix(i,:) = (filteredBinaryDischargeMatrix(i,:) ./ max(filteredBinaryDischargeMatrix(i,:))) .* max_IDR;
        end
    else
        error("Please select a valid normalization method")
    end
end

end

