function RUN_PHD()
% RUN_PHD runs the classic PHD filter. The following functions are called:
%
% PHD_initialiser   initialises the parameters of the filter
% simulator         might simulate a scenario if needed
% PHD_prediction    runs the PHD prediction
% PHD_update        runs the PHD update
%
% AUTHOR    Isabel Schlangen, (c) 2016

clear
clc
rng(5);

%% initialise constants and create Gaussian Mixture
cst = PHD_initialiser;

%% read folder where the code is located and create simulation folder
[homefolder,~,~] = fileparts(mfilename('fullpath'));
simfolder = [homefolder, '/simulation/'];

if ~isdir(simfolder)
    mkdir(simfolder);        %if the folder does not exist, create it
end

% [data,gt,trajectories] = simulator(simfolder,cst,0);
addpath '..'
[data,gt,trajectories] = dtSim(simfolder,cst,0);

%% create initial Gaussian structure
gaussian = struct('w', 0,...                    % weight
                  'm', zeros(cst.Nx,1),...      % mean
                  'C', eye(cst.Nx),...          % covariance 
                  'i', 0);                      % flag: active/inactive

gmm_u = repmat(gaussian,cst.gmmax,1);
isactive = zeros(1,cst.gmmax);

%% The main loop
for tt=1:cst.tmax
    [gmm_p,isactive] = PHD_prediction(gmm_u,isactive,cst);
    ind_p = find(isactive);
    [gmm_u,~,~,isactive] = PHD_update(gmm_p,data{tt},isactive,cst);
    ind_u = find(isactive);
    fprintf('time %3.d: #targets=%d, #meas=%d, pred - %3.d comp, mu=%.4g, update - %3.d comp, mu=%.4g \n',...
        tt,size(gt{tt},1),size(data{tt},1), length(ind_p),sum([gmm_p(ind_p).w]),length(ind_u),sum([gmm_u(ind_u).w]));
    plotGM(data{tt},gt{tt},trajectories,gmm_u,cst,tt);
end