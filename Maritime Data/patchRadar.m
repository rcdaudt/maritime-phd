function patchRadar
%{
    The radar detection file has detections as they come with per line scan
    of radar transmitter-recieiver cylce. This function converts the the
    detections to get a radar sweep to sync with the AIS data.
%}

load('tracksAIS_Run_01_005_patched.mat');
load('dataRadar/MarCE_Radar_Detections_Run_01_005.mat');

data = {};

for i = 1:numel(shipTime)-1
    idx = find(detections.time > shipTime(i)*1000  &  detections.time < shipTime(i+1)*1000);
    cData = struct();
    for j = 1:numel(idx)
        cData(j).t = detections.time(idx(j));
        cData(j).TR = [pi/2 - detections.azimuth(idx(j)); detections.range(idx(j))];
        cData(j).cov = [detections.covMatrices(idx(j),1) detections.covMatrices(idx(j),2); detections.covMatrices(idx(j),3) detections.covMatrices(idx(j),4)];
    end
    data = [data; cData];
end

save('MarCE_Radar_Detections_01_005_patched.mat','data')

% shipHeading = shipHeading(1:end-1);
shipTrajectoryX = shipTrajectoryX(1:end-1);
shipTrajectoryY = shipTrajectoryY(1:end-1);
shipTime = shipTime(1:end-1);


% save('tracksAIS_Run_01_005_patched.mat','shipHeading','shipTrajectoryX','shipTrajectoryY','shipTime');
save('tracksAIS_Run_01_005_patched.mat','shipTrajectoryX','shipTrajectoryY','shipTime');