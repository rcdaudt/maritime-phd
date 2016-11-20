function RUN_PHD_RealData_easy()
% RUN_PHD runs the classic PHD filter. The following functions are called:
%
% PHD_initialiser   initialises the parameters of the filter
% simulator         might simulate a scenario if needed
% PHD_prediction    runs the PHD prediction
% PHD_update        runs the PHD update
%
% AUTHOR    Isabel Schlangen, (c) 2016

clear
close all
clc
rng(5);

%% initialise constants and create Gaussian Mixture
cst = PHD_initialiser2;

%% read folder where the code is located and create simulation folder
[homefolder,~,~] = fileparts(mfilename('fullpath'));
simfolder = [homefolder, '/simulation/'];

if ~isdir(simfolder)
    mkdir(simfolder);        %if the folder does not exist, create it
end

addpath '..'
load('dtGroundTruthAIS.mat');
load('MarCE_Radar_Detections_01_005_patched_clean.mat');
%[data,gt,trajectories] = simulator(simfolder,cst,0);
% [data,gt,trajectories] = dtSim(simfolder,cst,0);

% ospa = Ospa_Adapted(data, gt, 1, 1);
% figure(); plot(ospa);
% title('Ospa metric');
%% create initial Gaussian structure
gaussian = struct('w', 0,...                    % weight
                  'm', zeros(cst.Nx,1),...      % mean
                  'C', eye(cst.Nx),...          % covariance 
                  'i', 0);                      % flag: active/inactive

gmm_u = repmat(gaussian,cst.gmmax,1);
isactive = zeros(1,cst.gmmax);

%% The main loop
% ospa = zeros(1,size(data,1));

offset = 1500; %2000;
scale = 1; %50/7500;
ospa = zeros(1,cst.tmax);
figure(420); hold on; box on; grid on; axis([0 4500 0 3000]);
gmm_u_s = [];

for tt=1:cst.tmax
    gt{tt} = (gt{tt}+offset)*scale;
    [gmm_p,isactive] = PHD_prediction(gmm_u,isactive,cst);
    ind_p = find(isactive);
    data = data_clean{tt};
    if numel(data) ~= 0
        TR_car = cell2mat({data.TR_car})';
        TR_car = (TR_car+offset)*scale;
        [gmm_u,~,~,isactive] = PHD_update(gmm_p,TR_car,isactive,cst);
    end
    
    ind_u = find(isactive);
    w = [gmm_u(ind_u).w];
    w_s = sort(w,'descend');
    n_obj = ceil(sum(w));
    ind_u_selected = ind_u(w >= w_s(min(n_obj,numel(w))));
    
    
    gmm_u_saved = gmm_u_s;
    gmm_u_s = gmm_u(ind_u_selected);
    
    % hack 2250
    aaa = [gmm_u_s.m];
    aaa = aaa(1,:);
    gmm_u_s = gmm_u_s(aaa ~= 2250);
    
    
    ospa(tt) = Ospa_Adapted(gmm_u_s, gt{tt}, 1, 1);
   
    fprintf('time %3.d: #targets=%d, #meas=%d, pred - %3.d comp, mu=%.4g, update - %3.d comp, mu=%.4g \n',...
        tt,size(gt{tt},1),size(TR_car,1), length(ind_p),sum([gmm_p(ind_p).w]),length(ind_u),sum([gmm_u(ind_u).w]));
    
%     plotGM2(TR_car,gt{tt},gmm_u_s,cst,tt);
    dunc_gmphd_plot(gmm_u_saved, gmm_u_s, 420, 2)
    pause(0.01)
end
%figure(); plot(ospa);title('OSPA metric for real data'); grid on;
% figure(); plot(ospa); title('Ospa metric'); grid on;