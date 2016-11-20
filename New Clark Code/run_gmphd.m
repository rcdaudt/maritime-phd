%% Clear and csstomize MATLAB 
clc;
clear;

set(0,'defaultfigurecolor',[1 1 1])
set(0,'DefaultFigureWindowStyle','docked');

figure(101); clf(101); axis([-500 500 -500 500]); box on; hold on;
figure(102); clf(102); axis([-500 500 -500 500]); box on; hold on;

%% GMPHD modelling params
model.birthRFS = .2;
model.falseAlarms.mean = 10;
model.prune_T = .1*(model.birthRFS/(5+model.falseAlarms.mean)); %;.1*power(10,-3);
model.merge_U = .5;
model.pD = .9;
model.pS = .95;
model.dT = 1;
model.noise_process = .3;
nSigma = 2;
model.noise_sensor = (3*nSigma)^2;
model.F = [1 model.dT;...
           0 1];
model.F = [model.F zeros(size(model.F));...
           zeros(size(model.F)) model.F];
model.Q = model.noise_process * [power(model.dT,4)/4 power(model.dT,3)/2;
                                 power(model.dT,3)/2 power(model.dT,2)];
model.Q = [ model.Q zeros(size(model.Q));
            zeros(size(model.Q)) model.Q];
model.H = [1 0 0 0;...
           0 0 1 0];
model.R = model.noise_sensor * [1^2 0;
                                0 1^2];
model.oSpaceVolume = 1000*1000;
model.falseAlarms.density = model.falseAlarms.mean/model.oSpaceVolume;

disp(model)

%% Generate simulated data
% try
%   load('measurements.mat');
% catch
%   generate_data();
%   load('measurements.mat');
% end
load('measurements_sigma3_lamda10.mat');

%% Plot simulated measurements
figure(101); box on; grid on; hold on;
title('Measurements');
nTime = numel(sensorMeasurements);
for j = 1:numel(sensorMeasurements)
    plot(sensorMeasurements{1,j}.xMeas,sensorMeasurements{1,j}.yMeas,'x','Color',[.2 .2 .2]);
end

%% plot groundtruth
figure(102); box on; grid on; hold on;
title('Tracks');
for j = 1:numel(groundTruth)
    h_102(1) = plot(groundTruth(j).track.x, groundTruth(j).track.y,'-k','LineWidth',1);
end
h_102(2) = plot(-500,500,'.k');
legend(h_102,'Groundtruth','GMPHD');

duration = numel(sensorMeasurements);
MTStateGT = cell(duration,1);
for ii = 1:duration
    for jj = 1:numel(groundTruth)
        MTStateGT{ii} = [MTStateGT{ii} [groundTruth(jj).track.x(ii) groundTruth(jj).track.y(ii)]'];
    end
end

%% structure for hypotheses and tracks
structHyp = struct(...
    'wk',-1,...                % Probability for the hypothesis to exist, keep -1 to lets functions know first iteration
    'mk',zeros(4,1),...        % Mean of the hypothesis
    'Pk',zeros(4),...          % Covariance of the hypothesis
    'Sk', zeros(4),...         % Innovation covariance
    'Kk', 0,...                % Kalman gain
    'neta', 0,...              % target obseravation H*mk              
    'tag', 0,...               % A zero tag means unassigned
    'entry', -1,...
    'status', 0);

HypP = structHyp;
ID = [];
Tracks.id = [];
Tracks.k = [];
Tracks.x = [];
nTracks = 0;

%% Filter
ospaCutOff = 200;
ospaArray = ones(1,duration)*ospaCutOff;

tic;
for k = 1:numel(sensorMeasurements)
    % Prediction
    HypP = gmphd_predict(HypP, model, sensorMeasurements{1,k},k);
    
    % Update
    HypP = gmphd_update( HypP, model, sensorMeasurements{1,k});
    
    % Prune and Merge
    HypP = gmphd_merge( HypP, model.prune_T, model.merge_U );
    
    % State Extraction
    MTState = []; % Multi target state
    wk = extractfield(HypP,'wk');
    disp(['sum of wk:' num2str(sum(wk)) ', nEstimated:' num2str(round(sum(wk)))])
    figure(102);hold on; box on; grid on;
    for i = 1:round(sum(wk))
        if(i>numel(wk))
            break;
        end
        plot(HypP(i).mk(1),HypP(i).mk(3),'.r','MarkerSize',10);
        MTState = [MTState [HypP(i).mk(1) HypP(i).mk(3)]'];
    end
    pause(.01);
    disp(['Iteration:' num2str(k)]);
    
    ospaArray(k) = ospa_dist(MTState, MTStateGT{k}, ospaCutOff, 1);
    figure(103); 
    plot(1:k, ospaArray(1:k), '-k', 'LineWidth', 2);
    xlabel('time');
    ylabel('OSPA dist.');
    drawnow;
    
    % Track Extraction
    cTracks = extractfield(HypP,'status');
    cTracks = HypP(find(cTracks==1));
    for i_cTracks = 1:numel(cTracks);
        if isempty(find(ID==cTracks(i_cTracks).tag))
            ID = [ID cTracks(i_cTracks).tag];
            nTracks = nTracks+1;
        end
    end
    if ~isempty(cTracks)
        cTracks_tag = extractfield(cTracks,'tag');
    end
    for i_Tracks = 1:nTracks
        idx = find(cTracks_tag==ID(i_Tracks));
        if isempty(idx)
            continue;
        end
%         cTracks(idx)
        Tracks(i_Tracks).id = ID(i_Tracks);
        Tracks(i_Tracks).k = [Tracks(i_Tracks).k k];
        Tracks(i_Tracks).x = [Tracks(i_Tracks).x cTracks(idx(1)).mk];
    end
end
runtime = toc;
disp('------------------');
disp(['Runtime: ' num2str(runtime) 'sec'])

figure(101);
for i = 1:numel(Tracks)
    plot(Tracks(i).x(1,:), Tracks(i).x(3,:), '-', 'Color', rand(1,3),'LineWidth',2)
end