% Simulator for use in RUN_PHD
% Rodrigo Daudt and Duncan Sommer

% Missing features (not necessarily necessary):
% Does not save data to 'pathName'
% Does not take bounds into account
% Does not remove objects when they go out of bounds
% Does not output trajectories for plotting

function [Y,X,tracks] = dtSim(pathName,cst,displ)
% Generate random state X and measurement Y for N time-steps
% X = [x dx/dt y dy/dt]
% Y = [x y]

% Assertions
N = cst.tmax;
assert(N == round(N))
assert(N > 3) % why not

% Parameters
dt = cst.T; % time-step
lambda = cst.nBirth; % poisson parameter for new objects
lambda_fa = cst.nFA; % poisson parameter for false alarms
P_survival = cst.pS; % probability of object survival
P_measure = cst.pD; % probability of measurement
% GET THIS FROM STRUCTURE
p_noise = 0.01; % position noise
assert(cst.sigmavel_x == cst.sigmavel_y)
v_noise = cst.sigmavel_x; % velocity noise
assert(cst.pixsize_x == cst.pixsize_y)
m_noise = cst.pixsize_x; % measurement noise

% Matrices
A = [1 dt 0 0;0 1 0 0;0 0 1 dt;0 0 0 1];
H = [1 0 0 0;0 0 1 0];

X = cell(N,1);
Y = cell(N,1);
tracks = cell(0,1);

n0 = round(4*rand()) + 2;

X{1} = [cst.xwidth*rand(1,n0);randn(1,n0);cst.ywidth*rand(1,n0);randn(1,n0)];
Y{1} = H*X{1};


for i = 2:N
    
    survival_var = rand(size(X{i-1},2),1);
    X_old = X{i-1};
    X_old_pruned = X_old(:,survival_var<P_survival);
    
    n = size(X_old_pruned,2);
    
    % Motion model
    X{i} = A*X_old_pruned + [p_noise*randn(1,n);v_noise*randn(1,n);p_noise*randn(1,n);v_noise*randn(1,n)];
    
    % Kill if out of bounds
    x = X{i}(1,:);
    y = X{i}(3,:);
    inds_keep = ((x < 0) + (x > cst.xwidth) + (y < 0) + (y > cst.ywidth)) == 0;
    X{i} = X{i}(:,inds_keep);
    
    % New objects
    nn = poissrnd(lambda);
    if nn > 0
        new_objects = [cst.xwidth*rand(1,nn);randn(1,nn);cst.ywidth*rand(1,nn);randn(1,nn)];
        X{i} = [X{i} new_objects];
    end
    
    Y_full = H*X{i};
    measurement_var = rand(size(X{i},2),1);
    Y{i} = Y_full(:,measurement_var<P_measure);
    Y{i} = Y{i} + m_noise*randn(size(Y{i}));
    
    % False alarms
    n_fa = poissrnd(lambda_fa);
    if n_fa > 0
        FA = [cst.xwidth*rand(1,n_fa); cst.ywidth*rand(1,n_fa)];
        Y{i} = [Y{i} FA];
    end
    
end

% Format the output matrices. 
% Get rid of velocity info and transpose each timeframe.
for j = 1:N
    x = X{j};
    X{j} = x([1 3],:)';
    Y{j} = Y{j}';
end

end