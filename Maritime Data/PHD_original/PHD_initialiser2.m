function cst = PHD_initialiser2()

% PHD_INITIALISER collects the allocation of all major constants of the
%                 filter and returns them as a single structure.
%
% OUT       cst             structure with all the constants
%
% AUTHOR    Isabel Schlangen, (c) 2016

%% dimensions of state and observation space
cst.Nx = 4;                         % dimension of state space
cst.Nz = 2;                         % dimension of observation space

%% flags and thresholds
cst.clut = 2;                       % clutter model: 1=NB, 2=Poisson
cst.dist = 2;                       % distance measure: 1=M, 2=H, 3=B

cst.gating = 1;                     % flag for gating
cst.tgate = chi2inv(0.99,cst.Nz);   % threshold for chi^2 gating, 99% conf.

cst.pruning = 1;                    % flag for pruning
cst.tprune = 1e-10;                  % threshold for pruning

cst.merging = 1;                    % flag for merging
if cst.dist==1                      % use Mahalanobis distance
    cst.tmerge = 4;
elseif cst.dist == 2                % use Hellinger distance
    cst.tmerge = 0.8;
elseif cst.dist ==3                 % use Bhattacharyya distance
    cst.merge = 1; %CHECK FOR GOOD VALUE!
end

%% Field of View
cst.xwidth = 50; %-1500                        % size of FoV in x direction (unit length)
cst.ywidth = 50; %-1200                        % size of FoV in y direction (unit length)
cst.npix_x = 50;                        % number of pixels in x direction
cst.npix_y = 50;                        % number of pixels in y direction
cst.pixsize_x = cst.xwidth/cst.npix_x;  % pixel size in x direction
cst.pixsize_y = cst.ywidth/cst.npix_y;  % pixel size in x direction
cst.npix = cst.npix_x*cst.npix_y;       % total number of pixels

%% Time constants
cst.tmax = 252;                     % number of frames
cst.T = 1;                          % time step length


%% all the other constants
cst.nFA = 1;                        % mean no of fa/frame
cst.pFA = cst.nFA/cst.npix;         % probability of false alarms (per px)

cst.pS = 0.999;                      % probability of survival
cst.pD = 0.95;                       % probability of detection

cst.sigmaz_x = cst.pixsize_x;       % measurement noise standard deviation
cst.sigmaz_y = cst.pixsize_y;       % measurement noise standard deviation

cst.sigmavel_x = 0.3;               % std deviation of initial velocity
cst.sigmavel_y = 0.3;               % std deviation of initial velocity

cst.q_x = 1;                      % accelleration noise
cst.q_y = 1;                      % accelleration noise

cst.gmmax = 500;                   % maximum number of Gaussians

%% Target birth
cst.nBirth = 1;                   % mean no of births / frame
cst.varBirth = 1.0000001; %cst.nBirth;          % in this case, the birth is Poisson
%pre-defined birth component: weight, mean, covariance, activity flag
cst.birthcomp = struct(...
    'w', cst.nBirth,...
    'm', [cst.xwidth/2;0;cst.ywidth/2;0],...
    'C', diag([cst.xwidth/2,cst.sigmavel_x,cst.ywidth/2,cst.sigmavel_y]).^2,...
    'i', 1);

%% Kalman initialisation
%Transition model: F and Q
cst.F = [1 cst.T 0 0; 0 1 0 0; 0 0 1 cst.T; 0 0 0 1];
cst.Q = [cst.q_x^2*0.25*cst.T^4, cst.q_x^2*0.5*cst.T^3, 0, 0;...
    cst.q_x^2*0.5*cst.T^3, cst.q_x^2*cst.T^2,  0, 0;...
    0, 0,  cst.q_y^2*0.25*cst.T^4, cst.q_y^2*0.5*cst.T^3;...
    0, 0,   cst.q_y^2*0.5*cst.T^3,     cst.q_y^2*cst.T^2;];
%Observation model: H and R
cst.H = [1 0 0 0; 0 0 1 0];
cst.R = diag([cst.sigmaz_x^2,cst.sigmaz_y^2]);

