function [dMU, d_per_time, tau, responsibleUnits]  = mudisplacement(psth, Fs, varargin)
P = inputParser;
addRequired(P, 'psth', @iscell)
addRequired(P, 'Fs', @isnumeric)
% addParameter(P, 'maxLag', 0.025, @(x) isa(x, 'struct'))
addParameter(P, 'maxLag', 0, @(x) isa(x, 'struct'))
parse(P, psth, Fs, varargin{:});

nSamp = cellfun(@(X) size(X,2), psth); % origin, for Demo
% nSamp = cellfun(@(X) size(X,1), psth);

% condition indices per sample
condIdx = cellfun(@(idx,T) idx*ones(1,T),...
    num2cell(1:length(psth))', num2cell(nSamp), 'uni',false);


% concatenate cells
psth = cat(2,psth{:});
condIdx = cat(2,condIdx{:});
nSampTot = size(psth,2);

% time lags
maxLagSamp = round(P.Results.maxLag*Fs);
lags = -maxLagSamp:maxLagSamp;

% anonymous functions
drPlus = @(R1,R2) max(0, max(R1-R2,[],1));
D = @(t) min([drPlus(psth(:,t),psth); drPlus(psth,psth(:,t))], [], 1);

dMU = zeros(1,nSampTot);
tau = zeros(2,nSampTot);
responsibleUnits = zeros(2,nSampTot);

progressBar = waitbar(0, "Displacement for current MU pair/population");
d_per_time = zeros(nSampTot,nSampTot);
for t = 1:nSampTot
    d_per_time(:,t) = D(t);
    if P.Results.maxLag == 0
        displacement = max(D(t));
    else
        [displacement, t2] = max(D(t));
        if displacement > 0
            
            % add lag to each time step
            times = {t, t2};
            tLagged = {t + lags, t2 + lags};
            for ii = 1:2
                tLagged{ii} = tLagged{ii}(tLagged{ii}>0 & tLagged{ii}<nSampTot);
                tLagged{ii} = tLagged{ii}(condIdx(tLagged{ii}) == condIdx(times{ii}));
            end
            
            % sample psths at each lagged time
            psthLagged = {psth(:,tLagged{1}), permute(psth(:,tLagged{2}),[1,3,2])};
            
            % optimize displacement over lags
            DLagged = min([drPlus(psthLagged{1}, psthLagged{2}); drPlus(psthLagged{2}, psthLagged{1})], [], 1);
            [~,tau1,tau2] = ind2sub(size(DLagged), find(DLagged==min(DLagged(:)),1));
            
            displacement = DLagged(1,tau1,tau2);
            tau(:,t) = [tau1; tau2];
            if displacement > 0
                R1 = psthLagged{1}(:,tau1);
                R2 = psthLagged{2}(:,tau2);
                [~,mu1] = max(R1-R2);
                [~,mu2] = max(R2-R1);
                responsibleUnits(:,t) = [mu1; mu2];
            end
        end
    end
    dMU(t) = displacement;
    waitbar(t/nSampTot, progressBar, "Displacement for current MU pair/population");
end
close(progressBar);
end