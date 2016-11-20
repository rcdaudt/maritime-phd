function model = genModel(pD, Lamda_fa, nSigma)
%GENMODEL Summary of this function goes here
%   Detailed explanation goes here
model.birthRFS = .1;
model.falseAlarms.mean = Lamda_fa;
model.prune_T = .1*(model.birthRFS/(5+model.falseAlarms.mean)); %;.1*power(10,-3);
model.merge_U = .5;
model.pD = pD;
model.pS = .95;
model.dT = 1;
model.noise_process = .2;
nSigma = nSigma;
model.noise_sensor = (2*nSigma)^2;
model.F = [1 model.dT;...
           0 1];
model.F = [model.F zeros(size(model.F));...
           zeros(size(model.F)) model.F];
model.Q = model.noise_process * [power(model.dT,4)/4 power(model.dT,3)/2;
                                 power(model.dT,3)/2 power(model.dT,2)];
model.Q = [ model.Q zeros(size(model.Q));
            zeros(size(model.Q)) model.Q];
model.H = [1 0 0 0;...
           0 0 1 0];
model.R = model.noise_sensor * [1^2 0;
                                0 1^2];
model.oSpaceVolume = 1000*1000;
% model.xRes = sqrt(model.R(1,1));
% model.yRes = sqrt(model.R(2,2));
% model.nCell = model.oSpace/(model.xRes*model.yRes)
model.falseAlarms.density = model.falseAlarms.mean/model.oSpaceVolume;

end

