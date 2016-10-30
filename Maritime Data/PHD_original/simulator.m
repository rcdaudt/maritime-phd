function [data,gt,allobjects] = simulator(pathName,cst,displ)
% SIMULATOR simulates a scenario using a nearly-constant velocity model.
%
% IN            pathName    name of path where simulation is saved
%               tmax        final time
%               displ       flag for display
%
% OUT           data        cell array containing simulated objects per frame
%               allobjects  cell array with ground truth per object
%
% AUTHOR        Isabel Schlangen, (c) 2016

% rng(24)
% tmax=1000;
% displ=1;

%% Parameters

%parameters
A = [cst.T^2/2, 0; cst.T, 0; 0, cst.T^2/2; 0, cst.T;];

%% Figure
if displ
    figure(1)
    clf
end

%% pre-allocate structures to store simulated data
objects = cell(0,1);            % object cell (per object)
allobjects = cell(0,1);         % all objects over time (for GT file)
data = cell(cst.tmax,1);        % data storage cell (per frame)
gt = cell(cst.tmax,1);          % gt storage cell (per frame)

%% Simulation
for tt=1:cst.tmax
    % deaths:
    indDeaths = [];
    for ii=1:length(objects)
        if rand>cst.pS
            indDeaths = [indDeaths,ii];
        end
    end
    allobjects = cat(1,allobjects,objects(indDeaths));
    objects(indDeaths) = [];
    
    % prediction of persisting targets (kill targets outside of FoV)
    indDeaths = [];
    for ii=1:length(objects)
        newstate = cst.F*objects{ii}(end,1:4)'+A*([cst.q_x^2;cst.q_y^2].*randn(2,1));
        if newstate(1)<0 || newstate(1)>cst.xwidth || newstate(3)<0 || newstate(3)>cst.ywidth
            indDeaths = [indDeaths,ii];
        else
            newstate = [newstate',tt];
            objects{ii} = [objects{ii};newstate];
        end
    end
    allobjects = cat(1,allobjects,objects(indDeaths));
    objects(indDeaths) = [];
    
   
    % create birth according to the chosen model
    if cst.nBirth==cst.varBirth %Poisson
        nBirths = poissrnd(cst.nBirth);
    elseif cst.nBirth<cst.varBirth %negative binomial
        betaval = cst.nBirth/(cst.varBirth-cst.nBirth);
        alphaval = betaval*cst.nBirth;
        nBirths = nbinrnd(alphaval, betaval/(1+betaval));
    else %binomial
        p = 1-cst.varBirth/cst.nBirth;
        n = cst.nBirth/p;
        nBirths = binornd(n,p);
    end
    
    for ii=1:nBirths
        if ~isempty(objects)
            objects = cat(1,objects,[rand*cst.xwidth, cst.sigmavel_x^2*randn,...
                rand*cst.ywidth, cst.sigmavel_y^2*randn,tt]);
        else
            objects = cell(1,1);
            objects{1} = [rand*cst.xwidth, cst.sigmavel_x^2*randn,...
                rand*cst.ywidth, cst.sigmavel_y^2*randn,tt];
        end
    end
    
    for ii = 1:length(objects)
        gt{tt} = [gt{tt}; [objects{ii}(end,1),objects{ii}(end,3)]];
        if rand<cst.pD
            data{tt} = [data{tt}; [objects{ii}(end,1),objects{ii}(end,3)] ...
                + [cst.sigmaz_x^2,cst.sigmaz_y^2].*randn(1,2)]; %is detected
        end
    end
       
    % clutter
    switch cst.clut
        case 1
            nclu = nbinrnd(r,p);
        case 2 
            nclu = poissrnd(cst.nFA);
    end    

    data{tt} = [data{tt};[cst.xwidth*rand(nclu,1),cst.ywidth*rand(nclu,1)]];
    
    %% Plot
    if displ
%         clf
        hold on
        axis equal
        for ii=1:length(objects)
            plot(objects{ii}(:,1),objects{ii}(:,3),'b');
        end
        if ~isempty(gt{tt})
            plot(gt{tt}(:,1),gt{tt}(:,2),'bo')
        end
        plot(data{tt}(:,1),data{tt}(:,2),'r.')
        drawnow
    end
end
if ~isempty(objects)
    allobjects = cat(1,allobjects,objects);
end

%% SAVE DATA
% GROUND TRUTH ____________________________________________________________
% Output file in xml format
fid = fopen([pathName,'gt.xml'],'wt');
fprintf(fid,'<?xml version="1.0" encoding="ISO-8859-1" standalone="no"?>\n');
fprintf(fid,'<root>\n');
t = datetime('now','TimeZone','local','Format','eee MMM dd hh:mm:ss Z y');
dateandtime = datestr(t,'ddd mmm dd HH:MM:SS yyyy',2000);
fprintf(fid,['<TrackContestISBI2012 generationDateTime="',...
    dateandtime,...
    '" info="http://bioimageanalysis.org/track/">\n']);
for ii = 1:length(allobjects)
    if size(allobjects{ii},1)>0
        fprintf(fid,'<particle>\n');
        for tt = size(allobjects{ii},1)
            fprintf(fid,'<detection t="%g" x="%g" y="%g" z="0.0"/>\n',...
                allobjects{ii}(tt,5),allobjects{ii}(tt,1),allobjects{ii}(tt,3));
        end
        fprintf(fid,'</particle>\n');
    end
end
fprintf(fid,'</TrackContestISBI2012>\n');
fprintf(fid,'</root>');
% mdf format
fid = fopen([pathName,'gt.mdf'], 'w');
fprintf(fid, 'MTrackJ 1.5.0 Data File\n');
fprintf(fid, 'Displaying true true true 1 1 0 0 100 10 0 0 0 1 1 0 0 true true true false\n');
fprintf(fid, 'Assembly 1 FF0000\n');
fprintf(fid, 'Cluster 1 FF0000\n');
for ii=1:length(allobjects)
    if size(allobjects{ii},1)>0
        fprintf(fid, 'Track %i FF0000 true\n',ii);
        for tt=1:size(allobjects{ii},1)
            fprintf(fid, 'Point %i %6.4f %6.4f %1.1f %1.1f %1.1f\n', tt,...
                allobjects{ii}(tt,1), allobjects{ii}(tt,3),...
                1.0, allobjects{ii}(tt,5), 1.0);
        end
    end
end
fprintf(fid, 'End of MTrackJ Data File\n');
fclose(fid);
% "measurements" __________________________________________________________
% xml format
fid = fopen([pathName,'measurements.xml'],'wt');
fprintf(fid,'<?xml version="1.0" encoding="ISO-8859-1" standalone="no"?>\n');
fprintf(fid,'<root>\n');
t = datetime('now','TimeZone','local','Format','eee MMM dd hh:mm:ss Z y');
dateandtime = datestr(t,'ddd mmm dd HH:MM:SS yyyy',2000);
fprintf(fid,['<TrackContestISBI2012 generationDateTime="',...
    dateandtime,...
    '" info="http://bioimageanalysis.org/track/">\n']);
for tt = 1:length(data)
    if size(data{tt},1)>0
        for ii = 1:size(data{tt},1)
            fprintf(fid,'<particle>\n');
            fprintf(fid,'<detection t="%g" x="%g" y="%g" z="0.0"/>\n',...
                tt,data{tt}(ii,1),data{tt}(ii,2));
            fprintf(fid,'</particle>\n');
        end
    end
end
fprintf(fid,'</TrackContestISBI2012>\n');
fprintf(fid,'</root>');

% mdf format
fid = fopen([pathName,'measurements.mdf'], 'w');
fprintf(fid, 'MTrackJ 1.5.0 Data File\n');
fprintf(fid, 'Displaying true true true 1 1 0 0 100 10 0 0 0 1 1 0 0 true true true false\n');
fprintf(fid, 'Assembly 1 FF0000\n');
fprintf(fid, 'Cluster 1 FF0000\n');
trackcount=0;
for ii=1:length(data)
    for tt=1:size(data{ii},1)
        trackcount=trackcount+1;
        fprintf(fid, 'Track %i FF0000 true\n',trackcount);
        fprintf(fid, 'Point 1 %6.4f %6.4f %1.1f %1.1f %1.1f\n',...
            data{ii}(tt,1), data{ii}(tt,2),...
            1.0, ii, 1.0);
    end
end
fprintf(fid, 'End of MTrackJ Data File\n');
fclose(fid);
end