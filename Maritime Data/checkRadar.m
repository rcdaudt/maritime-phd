clc;
clear;
close all
load('MarCE_Radar_Detections_01_005_patched.mat');
load('dtGroundTruthAIS.mat')



figure('Position',[300 700 700 650]);
for i = 1:numel(data)
    TR = extractfield(data{i},'TR');
    TR = reshape(TR,2,[]);
    Azimuth = TR(1,:);
    Range = TR(2,:);
    
    [X,Y] = pol2cart(Azimuth, Range);
    plot(X,Y,'.b');
    axis([-7000 7000 -7000 7000]);hold on;
    
    A = gt{i};
    scatter(A(:,1),A(:,2),'r.');
    
%     pause(.01);
end
grid on;
xlabel('X')
ylabel('Y')
title('Radar Data with AIS Ground Truth')

legend('Radar Data','AIS Data')