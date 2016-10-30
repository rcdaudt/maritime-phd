clc;
clear;
close all
load('MarCE_Radar_Detections_01_005_patched.mat');

for i = 1:numel(data)
    TR = extractfield(data{i},'TR');
    TR = reshape(TR,2,[]);
    Azimuth = TR(1,:);
    Range = TR(2,:);
    
    [X,Y] = pol2cart(Azimuth, Range);
    plot(X,Y,'.b');
%     axis([-7000 7000 -7000 7000]);hold on;
    axis([-10000 10000 -10000 10000]);hold on;
    pause(.01);
end