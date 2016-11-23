function RUN_PHD_RealData()
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
% rng(5);

%% initialise constants and create Gaussian Mixture
cst = PHD_initialiser3;

%% read folder where the code is located and create simulation folder
[homefolder,~,~] = fileparts(mfilename('fullpath'));
simfolder = [homefolder, '/simulation/'];

if ~isdir(simfolder)
    mkdir(simfolder);        %if the folder does not exist, create it
end

addpath '..'
load('dtGroundTruthAIS.mat');
load('MarCE_Radar_Detections_01_005_patched.mat');
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

offset = 7000;
scale = 1;
ospa = zeros(1,cst.tmax);

figure('Position',[300 700 700 650]);
figure(420); hold on; box on; grid on; axis([0 14000 0 14000]);
title('Data associations on GM-PHD');
gmm_u_s = [];

for tt=1:cst.tmax
    gt{tt} = (gt{tt}+offset)*scale;
    [gmm_p,isactive] = PHD_prediction(gmm_u,isactive,cst);
    ind_p = find(isactive);
    data_t = data{tt};
    TR_car = cell2mat({data_t.TR})';
    [X,Y] = pol2cart(TR_car(:,1),TR_car(:,2));
    TR_car = [X Y];
    TR_car = (TR_car+offset)*scale;
    [gmm_u,~,~,isactive] = PHD_update(gmm_p,TR_car,isactive,cst);
    
    ind_u = find(isactive);

    w = [gmm_u(ind_u).w];
    w_s = sort(w,'descend');
    n_obj = ceil(sum(w));
    ind_u_selected = ind_u(w >= w_s(min(n_obj,numel(w))));
    
    if exist('gmm_u_s', 'var')
        gmm_u_saved = gmm_u_s;
    end
    gmm_u_s = gmm_u(ind_u_selected);
    
    
    % hack 2
    aaa = [gmm_u_s.m];
    aaa = aaa(1,:);
    gmm_u_s = gmm_u_s(abs(aaa - 7000) > 0.1);
    
    
    
    
    ospa(tt) = Ospa_Adapted(gmm_u_s, gt{tt}, 450, 2);

    fprintf('time %3.d: #targets=%d, #meas=%d, pred - %3.d comp, mu=%.4g, update - %3.d comp, mu=%.4g \n',...
        tt,size(gt{tt},1),size(TR_car,1), length(ind_p),sum([gmm_p(ind_p).w]),length(ind_u),sum([gmm_u(ind_u).w]));
%     plotGM2(TR_car,gt{tt},gmm_u,cst,tt);

    
    figure(420); hold on; box on; grid on;
    if exist('gmm_u_saved', 'var')
        axis([0 14000 0 14000]);
        dunc_gmphd_plot(gmm_u_saved, gmm_u_s, 420, 2)
%         drawnow;
    end
end
% figure(); plot(ospa);

figure(); plot(ospa); title('Ospa metric for full data'); grid on;