function [ Hyp, MOL ] = gmphd_update( Hyp, model, sensorScan, MOL)
%GMPHD_UPDATE Summary of this function goes here
%   Detailed explanation goes here

nHyp = numel(Hyp);
wk = extractfield(Hyp ,'wk');
Dk_old = sum(wk);
MOL.T1 = exp(-model.pD*sum(wk));

% construction of PHD update components
for j = 1:nHyp
    Hyp(j).neta = model.H * Hyp(j).mk;
    Hyp(j).Sk = model.R + model.H * Hyp(j).Pk * model.H';
    Hyp(j).Kk = Hyp(j).Pk * model.H' * pinv(Hyp(j).Sk);
    Hyp(j).Pk = (eye(4) - Hyp(j).Kk * model.H) * Hyp(j).Pk;
end

% Array of PHD normaliztion terms
nofData = numel(sensorScan.xMeas);
tabNormalization = zeros(1,nofData);

% update
l_count = 0;
for i_obs = 1:numel(sensorScan.xMeas)
    l_count=l_count+1;
    z = [sensorScan.xMeas(i_obs);...
         sensorScan.yMeas(i_obs)];
    w_sum = 0;
    for j = 1:nHyp
        Hyp(l_count*nHyp+j).wk = model.pD*Hyp(j).wk *  mvnpdf(z, Hyp(j).neta,Hyp(j).Sk);
        % mvnpdf = exp(-.5*(z-Hyp(j).neta)' * pinv(Hyp(j).Sk) * (z-Hyp(j).neta))/sqrt(((2*pi)^numel(z))*det(Hyp(j).Sk))
        Hyp(l_count*nHyp+j).mk = Hyp(j).mk + Hyp(j).Kk*(z-Hyp(j).neta);
        Hyp(l_count*nHyp+j).Pk = Hyp(j).Pk;
        w_sum = w_sum + Hyp(l_count*nHyp+j).wk;
        Hyp(l_count*nHyp+j).tag = Hyp(j).tag;
        Hyp(l_count*nHyp+j).entry = Hyp(j).entry;
        Hyp(l_count*nHyp+j).status = Hyp(j).status;
    end
    for j = 1:nHyp
      Hyp(l_count*nHyp+j).wk = Hyp(l_count*nHyp+j).wk/(model.falseAlarms.density + w_sum);
    end
    tabNormalization(i_obs) = w_sum + model.falseAlarms.density;
end
% nHyp = l_count*nHyp + nHyp;
% 
% Missed Detections
sumWeightMiss = 0;
for j = 1:nHyp
    sumWeightMiss = sumWeightMiss + model.pD*Hyp(j).wk;
    Hyp(j).wk = (1 - model.pD) * Hyp(j).wk;
%     Hyp(j).mk = Hyp(j).mk;
%     Hyp(j).Pk = Hyp(j).Pk;
%     Hyp(j).tag = Hyp(j).tag;
end

MOL.T1 = exp(-sumWeightMiss);
MOL.T2 = sum(log(tabNormalization));
MOL.T = MOL.T1 * MOL.T2;

end

