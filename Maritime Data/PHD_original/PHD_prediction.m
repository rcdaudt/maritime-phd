function [gmm,isactive] = PHD_prediction(gmm,isactive,cst)

% PHD_PREDICTION executes the prediction step of a GM-PHD filter on a given
%                Gaussian Mixture.
%
% IN        gmmin       input mixture
%           cst         the constant structure
%
% OUT       gmmout      output mixture
%
% AUTHOR    Isabel Schlangen, (c) 2016

active = find(isactive);

%predict all active components
for ii=active
    gmm(ii).m = cst.F * gmm(ii).m;
    gmm(ii).C = cst.F * gmm(ii).C * cst.F' + cst.Q;
    gmm(ii).w = cst.pS * gmm(ii).w;
end

% initialise new-born component
inactive = find(~isactive,1,'first');
gmm(inactive) = cst.birthcomp; 
isactive(inactive) = 1;