function [dMNP, normBins, taus]  = mnpdispersion(psth, Fs, varargin)
P = inputParser;
addRequired(P, 'psth', @iscell)
addRequired(P, 'Fs', @isnumeric)
addParameter(P, 'pNorm', 1, @isnumeric)
addParameter(P, 'normEpsilon', 1, @isnumeric)
addParameter(P, 'maxLag', 0.025, @(x) isa(x, 'struct'))
addParameter(P, 'lagType', 'local', @(x) ischar(x) && ismember(x,{'local','global'}))
addParameter(P, 'deltaDrv', 25, @isnumeric)
parse(P, psth, Fs, varargin{:});

psth = P.Results.psth;

% motor pool norm
poolNorm = cellfun(@(X) sum(abs(X).^P.Results.pNorm, 1).^(1/P.Results.pNorm), psth, 'uni', false);

% psth derivative
if P.Results.maxLag > 0
    deltaDrv = P.Results.deltaDrv;
    tStart = 1+deltaDrv;
    psthDrv = cellfun(@(x) (x(:,tStart+deltaDrv:end) - x(:,tStart-deltaDrv:end-2*deltaDrv))/(2*deltaDrv), psth, 'uni',false);
    psthDrv = cellfun(@(x) [nan(size(x,1),tStart+deltaDrv-1),x], psthDrv, 'uni',false);
else
    psthDrv = psth;
end

% sample indices per condition
nSamp = cellfun(@length, poolNorm);
sampIdx = cellfun(@(T) 1:T, num2cell(nSamp), 'uni', false);

% condition indices per sample
condIdx = cellfun(@(idx,T) idx*ones(1,T),...
    num2cell(1:length(psth))', num2cell(nSamp), 'uni',false);

% concatenate cells
psth = cat(2,psth{:});
poolNorm = cat(2,poolNorm{:});
psthDrv = cat(2,psthDrv{:});
sampIdx = cat(2,sampIdx{:});
condIdx = cat(2,condIdx{:});

[nUnit,nSampTot] = size(psth);

% time lags
maxLagSamp = round(P.Results.maxLag*Fs);
lags = -maxLagSamp:maxLagSamp;

% bin pool values
edges = 0:P.Results.normEpsilon:(P.Results.normEpsilon+max(poolNorm));
[~,~,binIdx] = histcounts(poolNorm, edges);
nBins = length(edges)-1;

% results structure
normBins = zeros(1,nBins);
dMNP = zeros(1,nBins);
t = zeros(2,nBins);
if strcmp(P.Results.lagType,'local')
    taus = zeros(2,nBins,nUnit);
    tau = repmat({zeros(1,1,nUnit)},2,1);
else % global
    taus = zeros(2,nBins);
    tau = {0,0};
end

for iBin = 1:nBins
    
    % identify potential violations (times with similar norms)
    tQuery = find(binIdx==iBin);
    
    if ~isempty(tQuery)
        
        % difference in population state activity
        dr = zeros(length(tQuery));
        
        for ii = 1:length(tQuery)
            jj = (ii+1):length(tQuery);
            
            t1 = tQuery(ii);
            t2 = tQuery(jj);
            
            dr(ii,jj) = sum(abs(psth(:,t1) - psth(:,t2)).^P.Results.pNorm,1).^(1/P.Results.pNorm);
        end
        
        [maxi,maxj] = ind2sub(size(dr),find(dr==max(dr(:)),1));
        
        t1 = tQuery(maxi);
        t2 = tQuery(maxj);
        
        if P.Results.maxLag == 0
            dispersion = max(dr(:));
        else
            if strcmp(P.Results.lagType, 'local')
                
                % shifted states
                xCorrected = psth(:,[t1,t2]);
                tau = repmat({zeros(1,1,nUnit)},2,1);
                
                % try to move the state at t2 closer to the state at t1
                connectingVector = psth(:,t2) - psth(:,t1);
                temporalCorrection = -connectingVector ./ psthDrv(:,t2);  % This would be perfect if (1) the Taylor approximation were perfect and (2) we did not cap the temporal shift
                temporalCorrection(~isfinite(temporalCorrection)) = 0;
                
                % limit correction
                temporalCorrection = sign(temporalCorrection) .* min(round(abs(temporalCorrection)),maxLagSamp);
                temporalCorrection = max(temporalCorrection, find(condIdx==condIdx(t2),1)-t2);
                temporalCorrection = min(temporalCorrection, find(condIdx==condIdx(t2),1,'last')-t2);
                
                % Limit and apply the correction
                for ii = 1:size(psth,1)
                    xCorrected(ii,2) = psth(ii,t2+temporalCorrection(ii)); % go get actual data based on the delta in time approximated from the taylor series
                end
                
                % we will not keep the corrected state if it made things worse (which in theory could happen both because of the approximation and because we cap the temporal shift)
                if norm(xCorrected(:,2) - psth(:,t1),P.Results.pNorm) > norm(psth(:,t2) - psth(:,t1),P.Results.pNorm) % if we are further away than we were (in terms of the 1-norm)
                    xCorrected(:,2) = psth(:,t2);
                else
                    tau{2} = reshape(temporalCorrection,[1,1,nUnit]);
                end
                
                % try to move the state at t1 closer to the corrected state at t2
                connectingVectorNew = xCorrected(:,2) - psth(:,t1);
                temporalCorrection = connectingVectorNew ./ psthDrv(:,t1);
                temporalCorrection(~isfinite(temporalCorrection)) = 0;
                
                % limit correction
                temporalCorrection = sign(temporalCorrection).*min(round(abs(temporalCorrection)),maxLagSamp);
                temporalCorrection = max(temporalCorrection, find(condIdx==condIdx(t1),1)-t1);
                temporalCorrection = min(temporalCorrection, find(condIdx==condIdx(t1),1,'last')-t1);
                
                % Limit and apply the correction
                for ii = 1:size(psth,1)
                    xCorrected(ii,1) = psth(ii,t1+temporalCorrection(ii));  % go get actual data based on temporal shift that we estimated would be good
                end
                % As above, we will not keep the corrected state if it made things worse (which in theory could happen both because of the approximation and because we cap the temporal shift)
                if norm(xCorrected(:,2) - xCorrected(:,1),P.Results.pNorm) > norm(xCorrected(:,2) - psth(:,t1),P.Results.pNorm) % if we are further away than we were (in terms of the 1-norm)
                    xCorrected(:,1) = psth(:,t1);
                else
                    tau{1} = reshape(temporalCorrection,[1,1,nUnit]);
                end
                
                dispersion = norm(xCorrected(:,2) - xCorrected(:,1),P.Results.pNorm);
                
            else % global
                
                % bound lags
                tau = cell(2,1);
                for ii = 1:2
                    if ii == 1
                        tRef = t1;
                    else
                        tRef = t2;
                    end
                    tau{ii} = tRef + lags;
                    tau{ii} = tau{ii}(tau{ii}>=1 & tau{ii}<=nSampTot);
                    tau{ii} = tau{ii}(condIdx(tau{ii})==condIdx(tRef));
                    tau{ii} = tau{ii} - tRef;
                end
                
                % lagged firing rate differences
                drLag = abs(psth(:,t2+tau{2}) - permute(psth(:,t1+tau{1}),[1,3,2]));
                
                % pick best lag
                optLag = [0 0];
                
                % difference in population norm
                z = sum(drLag.^P.Results.pNorm,1).^(1/P.Results.pNorm);
                dispersion = min(z(:));
                [~,optLag(2),optLag(1)] = ind2sub(size(z),find(z==dispersion,1));
                
                % save time lags
                tau{1} = tau{1}(optLag(1));
                tau{2} = tau{2}(optLag(2));
                
            end
        end
        
        % save comparison values
        normBins(iBin) = edges(iBin);
        dMNP(iBin) = dispersion;
        t(:,iBin) = sampIdx([t1;t2]);
        taus(:,iBin,:) = cat(1,tau{:});
    end
end
end