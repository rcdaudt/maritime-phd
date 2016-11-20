function genSimData(nSigma, nLamda, measFileName)
%GENERATE_DATA Summary of this function goes here
%   Generate simulated data for GMPHD filter

figure(101); clf(101); hold on; axis([-500 500 -500 500]); box on; grid on;
figure(102); clf(102); hold on; axis([-500 500 -500 500]); box on; grid on;
figure(103); clf(103); hold on; axis([-500 500 -500 500]); box on; grid on;
figure(104); clf(104); hold on; axis([-500 500 -500 500]); box on; grid on;

% Initialize targets
% nSigma = 2;
% nLamda = 10;
nTargets = 3;
simTime = 100;

for i = 1:nTargets
    target(i).x = 0;
    target(i).y = 0;
    target(i).vx = 1;
    target(i).vy = 1;
    target(i).spawn = 'n';
    target(i).spawnTime = -1;
    target(i).spawnTarget.vx = 0;
    target(i).spawnTarget.vy = 0;
    target(i).death = -1;
end

% Initialize individual target properties
target(1).x = -100;
target(1).y = -400;
target(1).vx = 0;
target(1).vy = 7;
target(1).death = 0;

target(2).x = 100;
target(2).y = 400;
target(2).vx = 0;
target(2).vy = -7;
% target(2).spawn = 'y';
% target(2).spawnTime = 30;
% target(2).spawnTarget.vx = 2;
% target(2).spawnTarget.vy = -4;
target(2).x = 300;
target(2).y = 0;
target(2).vx = -7;
target(2).vy = 0;

target(3).x = -300;
target(3).y = -300;
target(3).vx = 6;
target(3).vy = 6;
target(3).spawn = 'n';
target(3).spawnTime = 60;%round(simTime/1.5);
target(3).spawnTarget.vx = -3;
target(3).spawnTarget.vy = 6;

% Initialize target states
for i = 1:nTargets
    target(i).state = [ target(i).x;...
                        target(i).vx;...
                        target(i).y;...
                        target(i).vy];
    groundTruth(i).track.x = [];
    groundTruth(i).track.y = [];
    groundTruth(i).track.t = [];
end

% Sampling time
dT = 1;

% State transition model
F = [1 dT;...
     0 1];
   
% Input control matrix or matrix for accelerative white gaussian noise
G = [dT^2; ...
    dT];
  
% For 2D systems
F = [F zeros(size(F)); zeros(size(F)) F];
G = [G zeros(size(G)); zeros(size(G)) G];

% noise factor for generating trajectories
pNoise = .001;

sensorMeasurements = {};
for i = 1:simTime
  sensor.xMeas = [];
  sensor.yMeas = [];
  sensor.amp = [];
  sensorImg = zeros(1000);
   
  for i_targetCounter = 1:nTargets
    figure(101); hold on;
    plot(target(i_targetCounter).state(1), target(i_targetCounter).state(3),'.b');
    groundTruth(i_targetCounter).track.x = [groundTruth(i_targetCounter).track.x    target(i_targetCounter).state(1)];
    groundTruth(i_targetCounter).track.y = [groundTruth(i_targetCounter).track.y    target(i_targetCounter).state(3)];
    groundTruth(i_targetCounter).track.t = [groundTruth(i_targetCounter).track.t    i];
    
    if(target(i_targetCounter).spawn=='y' & target(i_targetCounter).spawnTime==i)
      nTargets = nTargets + 1;
      target(nTargets).x = target(i_targetCounter).state(1);
      target(nTargets).y = target(i_targetCounter).state(3);
      target(nTargets).vx = target(i_targetCounter).spawnTarget.vx;
      target(nTargets).vy = target(i_targetCounter).spawnTarget.vy;
      target(nTargets).spawn = 'n';
      target(nTargets).spawnTime = -1;
      target(nTargets).spawnTarget.vx = 0;
      target(nTargets).spawnTarget.vy = 0;
      target(nTargets).state = [ target(nTargets).x;...
                          target(nTargets).vx;...
                          target(nTargets).y;...
                          target(nTargets).vy];
        groundTruth(nTargets).track.x = [];
        groundTruth(nTargets).track.y = [];
        groundTruth(nTargets).track.t = [];
    end
    
    SQRT = sqrt(target(i_targetCounter).state(2)^2 + target(i_targetCounter).state(4)^2);
    pNoiseX = pNoise * target(i_targetCounter).state(2)/SQRT + .05;
    pNoiseY = pNoise * target(i_targetCounter).state(4)/SQRT + .05;

    wk = [pNoiseX;pNoiseY] .* wgn(2,1,10);
    target(i_targetCounter).state = F * target(i_targetCounter).state+ G * wk;
    
    if rand < .9 % Pdetection
      sensor.xMeas = [sensor.xMeas target(i_targetCounter).state(1)+nSigma*randn];
      sensor.yMeas = [sensor.yMeas target(i_targetCounter).state(3)+nSigma*randn];
      sensor.amp = [sensor.amp round(55+200*rand)];
      imgIndex = sub2ind([1000 1000],501+round(sensor.yMeas(end)),1000-500+round(sensor.xMeas(end)));
      sensorImg(500+round(sensor.yMeas(end)),1000-500+round(sensor.xMeas(end))) = sensor.amp(end);
    end
  end
  
  % DEATH OF TARGETS
  killTargetIdx = [];
  for i_targetCounter = 1:nTargets
    if i==target(i_targetCounter).death
      killTargetIdx = [killTargetIdx i_targetCounter];
      nTargets = nTargets-1;
    end
  end
  target(killTargetIdx) = [];
  
%   figure(103); axis([0 1000 0 250]);
%   h_intensity_x(1) = plot(max(sensorImg,[],1),'Color',[0.9 0.5 0]);
%   figure(104); axis([0 1000 0 250]);
%   h_intensity_y(1) = plot(max(sensorImg,[],2),'Color',[0.9 0.5 0]);
%   sensorImg = zeros(1000);
  
  PRN = poissrnd(nLamda);
%   disp(['PRN is ' num2str(PRN)]);
  for i_randCounter = 1:PRN
    sensor.xMeas = [sensor.xMeas -500+1000*rand];
    sensor.yMeas = [sensor.yMeas -500+1000*rand];
    sensor.amp = [sensor.amp round(127*rand)];
    %imgIndex = sub2ind([1000 1000],500+round(sensor.yMeas(end)),1000-500+round(sensor.xMeas(end)));
    sensorImg(501+round(sensor.yMeas(end)),1001-500+round(sensor.xMeas(end))) = sensor.amp(end);
  end
  % Add random noise points
%   for i_randCounter = 1:floor(rand*5)
%     sensor.xMeas = [sensor.xMeas -500+1000*rand];
%     sensor.yMeas = [sensor.yMeas -500+1000*rand];
%   end
  sensorMeasurements{i} = sensor;
  
  figure(101); hTitle_101 = title('Target trajectories');
  figure(102); hold on; hTitle_102 = title('Sensor observation');
  plot(sensor.xMeas,sensor.yMeas,'.r');
%   figure(103); hold on; hTitle_x = title('Intensity in x-axis'); ylabel('intensity');axis([0 1000 0 250]);
%   h_intensity_x(2) = plot(max(sensorImg,[],1),'Color',[0 0.5 .8]); hold off;
%   legend(h_intensity_x,'Target','Noise');
%   figure(104); hold on; hTitle_y = title('Intensity in y-axis'); ylabel('intensity');axis([0 1000 0 250]);
%   h_intensity_y(2) = plot(max(sensorImg,[],2),'Color',[0 0.5 .8]); hold off;
%   legend(h_intensity_y,'Target','Noise');
  
  figure(101);
  % Adjust font and axes properties
    set( gca                       , ...
        'FontName'   , 'Helvetica' );
    set( hTitle_101                    , ...
        'FontSize'   , 10          , ...
        'FontWeight' , 'bold'      );
  figure(102);
  % Adjust font and axes properties
    set( gca                       , ...
        'FontName'   , 'Helvetica' );
    set( hTitle_102                    , ...
        'FontSize'   , 10          , ...
        'FontWeight' , 'bold'      );
%   figure(103);
%   % Adjust font and axes properties
%     set( gca                       , ...
%         'FontName'   , 'Helvetica' );
%     set( hTitle_x                    , ...
%         'FontSize'   , 10          , ...
%         'FontWeight' , 'bold'      );
%   figure(104);
%   % Adjust font and axes properties
%     set( gca                       , ...
%         'FontName'   , 'Helvetica' );
%     set( hTitle_y                    , ...
%         'FontSize'   , 10          , ...
%         'FontWeight' , 'bold'      );
% 
%   if(i==1)
%     figure(101);
%     outfile = 't1.gif';
%     frame = getframe(101);
%     frame = frame.cdata;
%     [imind,cm] = rgb2ind(frame,256);
%     imwrite(imind,cm,outfile,'gif','DelayTime',0.05,'loopcount',1);
%         figure(102);
%     outfile = 't2.gif';
%     frame = getframe(102);
%     frame = frame.cdata;
%     [imind,cm] = rgb2ind(frame,256);
%     imwrite(imind,cm,outfile,'gif','DelayTime',0.05,'loopcount',1);
%         figure(103);
%     outfile = 't3.gif';
%     frame = getframe(103);
%     frame = frame.cdata;
%     [imind,cm] = rgb2ind(frame,256);
%     imwrite(imind,cm,outfile,'gif','DelayTime',0.05,'loopcount',1);
%         figure(104);
%     outfile = 't4.gif';
%     frame = getframe(104);
%     frame = frame.cdata;
%     [imind,cm] = rgb2ind(frame,256);
%     imwrite(imind,cm,outfile,'gif','DelayTime',0.05,'loopcount',1);
%   else
%     figure(101)
%     frame = getframe(101);
%     frame = frame.cdata;
%     [imind,cm] = rgb2ind(frame,256);
%     outfile = 't1.gif';
%     imwrite(imind,cm,outfile,'gif','DelayTime',0.05,'writemode','append');
%         figure(102)
%     frame = getframe(102);
%     frame = frame.cdata;
%     [imind,cm] = rgb2ind(frame,256);
%     outfile = 't2.gif';
%     imwrite(imind,cm,outfile,'gif','DelayTime',0.05,'writemode','append');
%         figure(103)
%     frame = getframe(103);
%     frame = frame.cdata;
%     [imind,cm] = rgb2ind(frame,256);
%     outfile = 't3.gif';
%     imwrite(imind,cm,outfile,'gif','DelayTime',0.05,'writemode','append');
%         figure(104);
%     frame = getframe(104);
%     frame = frame.cdata;
%     [imind,cm] = rgb2ind(frame,256);
%     outfile = 't4.gif';
%     imwrite(imind,cm,outfile,'gif','DelayTime',0.05,'writemode','append');
%   end5
  
  pause(.01);
  sensorImg = zeros(1000);
end

save( measFileName, 'sensorMeasurements', 'groundTruth');

% disp('Press any key to continue'); pause;
% figure(101); clf(101); box on; grid on; hold on;
% figure(102); clf(102); box on; grid on; hold on;

end

