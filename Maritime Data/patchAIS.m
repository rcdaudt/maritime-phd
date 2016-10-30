function patchAIS
%{
    This function loads the MarCE AIS groundtruth data and extracts just
    the data of own ship. It syncs this data in accordance to the sampling
    time of each radar sweep. 
%}

load('dataAIS/tracksAIS_Run_01_005.mat');
ownTrack  = tracksAIS(1); clear('tracksAIS');
startTime = ownTrack.numTime(1);
endTime = ownTrack.numTime(end);

% sampling rate, time it takes to get one radar sweep
model.T = 2.4;                        

ownData.T = startTime:model.T:endTime;
ownData.X = zeros(1,numel(ownData.T));
ownData.Y = zeros(1,numel(ownData.T));
ownData.H = zeros(1,numel(ownData.T));
ownData.startTime = startTime;
ownData.endTime = endTime;
duration = numel(ownData.T);

for i = 1:duration
    [~,idx] = min(abs(ownTrack.numTime-ownData.T(i)));
    ownData.X(i) = ownTrack.X(idx);
    ownData.Y(i) = ownTrack.Y(idx);
    ownData.H(i) = ownTrack.Heading(idx);
end

% shipHeading = ownData.H;
shipTrajectoryX = ownData.X;
shipTrajectoryY = ownData.Y;
shipTime = ownData.T;


% save('tracksAIS_Run_01_005_patched.mat','shipHeading','shipTrajectoryX','shipTrajectoryY','shipTime');
save('tracksAIS_Run_01_005_patched.mat','shipTrajectoryX','shipTrajectoryY','shipTime');