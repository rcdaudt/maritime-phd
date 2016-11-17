function dtExtractGroundTruth


load('tracksAIS_Run_01_005_patched.mat');
load('dataAIS/tracksAIS_Run_01_005.mat');

ship1 = tracksAIS(5);
ship2 = tracksAIS(25);


gt = cell(numel(shipTime),1);

for i = 1:numel(shipTime)
    idx1 = find(ship1.numTime > shipTime(i), 1);
    idx2 = find(ship2.numTime > shipTime(i), 1);
    
    coords = [ship1.X(idx1) ship2.X(idx2);...
              ship1.Y(idx1) ship2.Y(idx2)]';
    
    
    
    gt{i} = coords;
end


save('dtGroundTruthAIS.mat','gt')


% save('tracksAIS_Run_01_005_patched.mat','shipHeading','shipTrajectoryX','shipTrajectoryY','shipTime');
% save('tracksAIS_Run_01_005_patched.mat','shipTrajectoryX','shipTrajectoryY','shipTime');