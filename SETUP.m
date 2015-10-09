function Params = SETUP
%% Load library
run('./../SET_PATH.m')
set_dpopt;
% run CMEX file to preload memory
pp = struct('form','MKLpp','breaks',{{[1 2 3 4]}},...
    'Values',[1 2 3 4],'coefs',[],'order',[4],...
    'Method',[],'ExtrapolationOrder',[],'thread',1,...
    'orient','curvefit');
warmUpPp = myppual(pp);

%% Parameters
% utility
Sigma = 2;
Chi = 0.3;
Beta = 0.98057846155234;
USigma = 0.3;

% production
Gamma = 0.36;
Theta = 0.71;
Delta = 0.0157;

% Tax
TauPiBar = 0;
TauRBar = 0;
TauLBar = 0;

%
ZBar = 1;
RBar = 0.05/4;

% Shock
% Aggregate
ZBar = 1;
ZRho = 0.75;
ZSigma = 0.012;
% G
GRho = 0.8909;
GBar = 0.20;
BBar = 0.3597*4;
Rho1 = 0.017;
Rho2 = 0.484;
% Lambda
LambdaBar = 2.14;
LambdaRho = 0.9170;
LambdaSigma = 0.049;

% Epsilon
EpsilonPts = 3;
SuperEpsilon = 0;
EpsilonRho = 0.755^0.25;
EpsilonSigma = 0.1;

% 
ZetaPts = 7;
ZetaMu = 0.4172;
ZetaSigma = 0.3843;
ZetaRho = 0.7;
% Entre shock
EntreShock1 = 0.0023;
EntreShock2 = 0.9808;

% computation
NumOfThreads = 8;
TolOpt = 1e-14;
TolVfi = 1e-12;
TolEqSs = 1e-12;
ShowDetail = 0;

% Beta = 0.979679791267919;
CaliX0 = [0.980137707974586     0.020513      0.30921    0.0039207      0.97187      0.75967 ...
    0.22786      0.76197       1.3213       1.6862      0.27075];
x = CaliX0;
Beta = x(1);
Delta = x(2);
Chi = x(3);
EntreShock1 = x(4);
EntreShock2 = x(5);
ZetaMu = x(6);
ZetaSigma = x(7);
Theta = x(8);
LambdaBar = x(9);
BBar = x(10);
GBar = x(11);
clear x;

GSigma = ((1-GRho^2)*0.0010243)^0.5*GBar;

% Beta = 0.981610878232992;

Params = v2struct;
end