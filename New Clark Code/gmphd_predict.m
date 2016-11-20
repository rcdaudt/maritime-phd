function Hyp = gmphd_predict( Hyp, model, sensorScan,k )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nHyp = numel(Hyp);
if Hyp(1).wk==-1
    nHyp = 0;
end

% Prediction of existing targets
for j = 1:nHyp
    Hyp(j).wk = model.pS * Hyp(j).wk;
    Hyp(j).mk = model.F * Hyp(j).mk;
    Hyp(j).Pk = model.Q + model.F * Hyp(j).Pk * model.F';
end

tags = extractfield(Hyp,'tag');
nTags = max(tags);

%% Prediction of new births
% [x_pos,y_pos] = meshgrid([-250 250]);
% for j = 1:4
%     Hyp(nHyp+j).wk = .5 + abs(randn)/1000;
%     Hyp(nHyp+j).mk = [x_pos(j);...
%                        5;...
%                        y_pos(j);...
%                        5];
%     Hyp(nHyp+j).Pk = 27000* eye(4).*diag([1 .01 1 .01]);
% %     figure(102); hold on;
% %     h_ellips(j) = ellips(x_pos(j),y_pos(j),diag([Hyp(nHyp+j).Pk(1,1) Hyp(nHyp+j).Pk(3,3)]),'r');
% %     h_ellips(j) = ellips(0,0,diag([20000 20000]),'r');
% end

%% Measurement driven new births
for j = 1:numel(sensorScan.xMeas)
    Hyp(nHyp+j).wk = model.birthRFS/numel(sensorScan.xMeas) + abs(randn)/1000;
    Hyp(nHyp+j).mk = [sensorScan.xMeas(j);...
                      0;...
                      sensorScan.yMeas(j);...
                      0];
    Hyp(nHyp+j).Pk = 27000*eye(4).*diag([1 .001 1 .001]);
    % Giving a new tag to new target births
    Hyp(nHyp+j).entry = k;
    Hyp(nHyp+j).status = 0;
    Hyp(nHyp+j).tag = nTags+1;
    nTags = nTags + 1;
end

%% Measurement driven birth only at start
% if k == 1
% for j = 1:numel(sensorScan.xMeas)
%     Hyp(nHyp+j).wk = .1 + abs(randn)/1000;
%     Hyp(nHyp+j).mk = [sensorScan.xMeas(j);...
%                       0;...
%                       sensorScan.yMeas(j);...
%                       7];
%     Hyp(nHyp+j).Pk = eye(4).*diag([1 .001 1 .001]);
%     if j==2
%         Hyp(nHyp+j).mk(4) = -7;
%     end
% end
% end

end

