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
ZRho = 0.9;
ZSigma = 0.01;
% G
GRho = 0.89;
GSigma = 0.02;
GBar = 0.20;
BBar = 0.3597*4;
Rho1 = 0.017;
Rho2 = 0.484;
% Lambda
LambdaBar = 2.14;
LambdaRho = 0.9;
LambdaSigma = 0.1;

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
TolOpt = 1e-12;
TolVfi = 1e-8;
TolEqSs = 1e-8;
ShowDetail = 0;

CaliX0 = [0.980795682546104    0.0223    0.2914    0.0025    0.9760    0.2388    0.5826    0.3264 ...
    0.7600    1.2119    1.6424    0.2614];
x = CaliX0;
Beta = x(1);
Delta = x(2);
Chi = x(3);
EntreShock1 = x(4);
EntreShock2 = x(5);
USigma = x(6);
ZetaMu = x(7);
ZetaSigma = x(8);
Theta = x(9);
LambdaBar = x(10);
BBar = x(11);
GBar = x(12);
clear x;

Beta = 0.981610878232992;

Params = v2struct;
end