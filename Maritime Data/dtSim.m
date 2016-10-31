function [X,Y] = dtSim(N)
% Generate random state X and measurement Y for N time-steps
% X = [x dx/dt y dy/dt]
% Y = [x y]
% Rodrigo Daudt

% Assertions
assert(N == round(N))
assert(N > 3) % why not

% Parameters
dt = 0.1; % time-step
lambda = 0.3; % poisson parameter for new objects
lambda_fa = 1; % poisson parameter for false alarms
p_noise = 0.01; % position noise
v_noise = 0.1; % velocity noise
m_noise = 0.05; % measurement noise
P_survival = 0.97; % probability of object death
P_measure = 0.9; % probability of measurement

% Matrices
A = [1 dt 0 0;0 1 0 0;0 0 1 dt;0 0 0 1];
H = [1 0 0 0;0 0 1 0];

X = cell(N,1);
Y = cell(N,1);


n0 = round(4*rand()) + 2;

X{1} = [rand(1,n0);randn(1,n0);rand(1,n0);randn(1,n0)];
Y{1} = H*X{1};


for i = 2:N
    
    survival_var = rand(size(X{i-1},2),1);
    X_old = X{i-1};
    X_old_pruned = X_old(:,survival_var<P_survival);
    
    n = size(X_old_pruned,2);
    
    % Motion model
    X{i} = A*X_old_pruned + [p_noise*randn(1,n);v_noise*randn(1,n);p_noise*randn(1,n);v_noise*randn(1,n)];
    
    % New objects
    nn = poissrnd(lambda);
    if nn > 0
        new_objects = [rand(1,nn);randn(1,nn);rand(1,nn);randn(1,nn)];
        X{i} = [X{i} new_objects];
    end
    
    Y_full = H*X{i};
    measurement_var = rand(size(X{i},2),1);
    Y{i} = Y_full(:,measurement_var<P_measure);
    Y{i} = Y{i} + m_noise*randn(size(Y{i}));
    
    % False alarms
    n_fa = poissrnd(lambda_fa);
    if n_fa > 0
        Y{i} = [Y{i} rand(2,n_fa)];
    end
    
end



end