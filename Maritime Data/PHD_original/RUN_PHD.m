function RUN_PHD()
% RUN_PHD runs the classic PHD filter. The following functions are called:
%
% PHD_initialiser   initialises the parameters of the filter
% simulator         might simulate a scenario if needed
% PHD_prediction    runs the PHD prediction
% PHD_update        runs the PHD update
%
% AUTHOR    Isabel Schlangen, (c) 2016



clear all
close all
clc

% rng(444);
rng(1);

dt_tic = tic;

nruns = 256;
ospa_acc = zeros(1,200);

%% initialise constants and create Gaussian Mixture
cst = PHD_initialiser4;

for i_run = 1:nruns

    %% read folder where the code is located and create simulation folder
    % [homefolder,~,~] = fileparts(mfilename('fullpath'));
    % simfolder = [homefolder, '/simulation/'];
    % 
    % if ~isdir(simfolder)
    %     mkdir(simfolder);        %if the folder does not exist, create it
    % end

    addpath '..'
    [data,gt,trajectories] = dtSim('holaquetal',cst,0);
    % dtWatevs3(data)
    % dtWatevs4(trajectories)
    % drawnow;


    %% create initial Gaussian structure
    gaussian = struct('w', 0,...                    % weight
                      'm', zeros(cst.Nx,1),...      % mean
                      'C', eye(cst.Nx),...          % covariance 
                      'i', 0);                      % flag: active/inactive

    gmm_u = repmat(gaussian,cst.gmmax,1);
    isactive = zeros(1,cst.gmmax);

    %% The main loop
    ospa = zeros(1,size(data,1));

%     figure(420); hold on; box on; grid on; axis([0 50 0 50]);
%     title('Data associations on GM-PHD');
    % title('Estimated targets on GM-PHD');

    for tt=1:cst.tmax
    % for tt=1:25
        [gmm_p,isactive] = PHD_prediction(gmm_u,isactive,cst);
        ind_p = find(isactive);
        [gmm_u,~,~,isactive] = PHD_update(gmm_p,data{tt},isactive,cst);

        ind_u = find(isactive);

        w = [gmm_u(ind_u).w];
        w_s = sort(w,'descend');
        n_obj = ceil(sum(w));
        if numel(w) == 0
            continue;
        end
        ind_u_selected = ind_u(w >= w_s(min(n_obj,numel(w))));


        if exist('gmm_u_s', 'var')
            gmm_u_saved = gmm_u_s;
        end
        gmm_u_s = gmm_u(ind_u_selected);


        % hack 2
        aaa = [gmm_u_s.m];
        aaa = aaa(1,:);
        gmm_u_s = gmm_u_s(aaa ~= 25);


        ospa(tt) = Ospa_Adapted(gmm_u_s, gt{tt}, 1, 2);


        fprintf('time %3.d: #targets=%d, #meas=%d, pred - %3.d comp, mu=%.4g, update - %3.d comp, mu=%.4g \n',...
            tt,size(gt{tt},1),size(data{tt},1), length(ind_p),sum([gmm_p(ind_p).w]),length(ind_u),sum([gmm_u(ind_u).w]));
    %     plotGM(data{tt},gt{tt},trajectories,gmm_u,cst,tt);

%         figure(420); hold on; box on; grid on;
%         if exist('gmm_u_saved', 'var')
%             axis([0 50 0 50]);
%             dunc_gmphd_plot(gmm_u_saved, gmm_u_s, 420, 2)
%     %         drawnow;
%         end
    end
    
    ospa_acc = ospa_acc + ospa;
end


figure(); plot(ospa_acc/nruns); title('Ospa metric'); grid on;

toc(dt_tic)