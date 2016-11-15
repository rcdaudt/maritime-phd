function [] = dtFindBoats()

% Good boats seem to be 5 and 25

clear all
close all
clc


load 'dataAIS/tracksAIS_Run_01_005.mat'

figure;

for i = 1:numel(tracksAIS)
    X = tracksAIS(i).X;
    Y = tracksAIS(i).Y;

    plot(X,Y);
    grid on;
    
    pause(0.01);
end


end