% Trim the good ship data from the complete dataset
% Rodrigo Daudt

clc;
clear;
load('MarCE_Radar_Detections_01_005_patched.mat');
load('dtGroundTruthAIS.mat')

data_clean = cell(size(data));

for i = 1:numel(data)-10
    TR = extractfield(data{i},'TR');
    TR = reshape(TR,2,[]);
    Azimuth = TR(1,:);
    Range = TR(2,:);
    
    [X,Y] = pol2cart(Azimuth, Range);
    
    ind = 1;
    for j = 1:length(X)
        if (X(j) > -1500) && (X(j) < 2800) && (Y(j) > -1200) && (Y(j) < 1500)
            data_clean{i}(ind).t = data{i}(j).t;
            data_clean{i}(ind).TR = data{i}(j).TR;
            data_clean{i}(ind).TR_car = [X(j);Y(j)];
            data_clean{i}(ind).cov = data{i}(j).cov;
            ind = ind + 1;
        end
    end
end

figure('Position',[300 700 700 650]);
for i = 1:numel(data_clean)
    try
        TR = extractfield(data_clean{i},'TR');
    end
    TR = reshape(TR,2,[]);
    Azimuth = TR(1,:);
    Range = TR(2,:);
    
    [X,Y] = pol2cart(Azimuth, Range);
    plot(X,Y,'.b');
    grid on;
    axis([-1500 2800 -1200 1500]);hold on;
    
    A = gt{i};
    scatter(A(:,1),A(:,2),'r.');
%     axis([-10000 10000 -10000 10000]);hold on;
%     pause(.01);
end
grid on;
xlabel('X')
ylabel('Y')
title('Clean Radar Data with AIS Groundtruth')

legend('Radar Data','AIS Data')

% save('MarCE_Radar_Detections_01_005_patched_clean.mat','data_clean')
