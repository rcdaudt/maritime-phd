function [gmm,variance,likelihood,isactive] = PHD_update(gmm,z,isactive,cst)

% PHD_UPDATE executes the update step of a GM-PHD filter on a given
%            Gaussian Mixture.
%
% IN        gmm       input mixture
%           cst       the constant structure
%
% OUT       gmm       updated mixture
%
% AUTHOR    Isabel Schlangen, (c) 2016

Nmeas = size(z,1);                % number of measurements
active_pred = find(isactive);     % collect all components

gmmtable = repmat(gmm(active_pred),1,Nmeas);

%% update all active components
normtable = cst.pFA*ones(Nmeas,1);
% ASSOCIATIONS
for jj=1:Nmeas
    for ii=1:length(active_pred)
        % innovation:
        y = z(jj,:)'-cst.H*gmm(active_pred(ii)).m;
        S = cst.H * gmm(active_pred(ii)).C * cst.H' + cst.R;

        % Kalman gain:
        K = (gmm(active_pred(ii)).C * cst.H')/S;
        
        gate = y'/S*y;
        
        % discard the combination if we want gating and m is outside of the
        % threshold around z:
        if ~(cst.gating && gate>cst.tgate) 
            if cst.Nz == 2
                detS = S(1,1)*S(2,2)-S(2,1)*S(1,2);
            else
                detS = det(S);
            end        
            q = exp(-0.5 * gate)/sqrt((2*pi)^(cst.Nz)*abs(detS));

            %update the Gaussian:
            gmmtable(ii,jj).w = q * cst.pD * gmm(active_pred(ii)).w;
            gmmtable(ii,jj).m = gmm(active_pred(ii)).m + K*y;
            gmmtable(ii,jj).C = (eye(cst.Nx)-K*cst.H) * gmm(active_pred(ii)).C;
            gmmtable(ii,jj).i = 1;
        else
            gmmtable(ii,jj).w = 0;
            gmmtable(ii,jj).i = 0;
        end
    end
    % normalise per measurement
    normtable(jj) = cst.pFA + sum([gmmtable(:,jj).w]);
    for ii=1:length(active_pred)
        gmmtable(ii,jj).w = gmmtable(ii,jj).w/normtable(jj);
        if gmmtable(ii,jj).w > 1
            keyboard;
        end
    end
end
if sum(sum([gmmtable(:).w]))>Nmeas
    keyboard;
end

% MISSED DETECTIONS
mol_mdterm = 0;
for ii=active_pred
    mol_mdterm = mol_mdterm + cst.pD*gmm(ii).w; %for the likelihood
    gmm(ii).w = (1-cst.pD)*gmm(ii).w;           %for weight update
end

likelihood = exp(-mol_mdterm)*prod(normtable);

%% variance
allweights =  zeros(size(gmmtable));
for ii=1:size(gmmtable,1)
   allweights(ii,:) = [gmmtable(ii,:).w]; 
end
variance = sum([gmm(ii).w]) + sum(sum(allweights.*(1-allweights)));

%% save all components that are still active for the next iteration
inactive = find(~isactive);  % collect indices of free space in gmm
activateind = 0;
for ii=1:length(active_pred)
    for jj=1:Nmeas
        if gmmtable(ii,jj).w > cst.pruning*cst.tprune
            activateind = activateind+1;
            if activateind>length(inactive)
                error('Maximum number of Gauss components reached.\n');            
            end
            gmm(inactive(activateind)) = gmmtable(ii,jj);
            isactive(inactive(activateind)) = 1;
        end
    end
end

%% Pruning
active_up = find(isactive);
prunedweights = 0;
for ii=active_up
    if gmm(ii).w<cst.pruning*cst.tprune
        prunedweights = prunedweights + gmm(ii).w;
        gmm(ii).w = 0;
        gmm(ii).i = 0;
        isactive(ii) = 0;
    end
end
active_pruned = find(isactive);
ncomp = length(active_pruned);
for ii=active_pruned
    gmm(ii).w = gmm(ii).w + prunedweights/ncomp;
end

% fprintf('MD: %g AS: %g, pruned:%g, total:%g\n',sum([gmm(active_pred).w]),sum([gmm(setdiff(active_up,active_pred)).w]),prunedweights,sum([gmm([gmm.i]==1).w]));

%% Merging
[gmm,isactive] = merging(gmm,isactive,cst);

% fprintf('total weight after merging: %g\n', sum([gmm([gmm.i]==1).w]));

