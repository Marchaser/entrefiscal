function NewParams = COMMON(Params)
v2struct(Params);
clear Params;

[EpsilonTrans,EpsilonGrid] = markovappr(EpsilonRho,EpsilonSigma,1,EpsilonPts);
EpsilonGrid = exp(EpsilonGrid);

%{
[ZetaTrans,ZetaGrid] = markovappr(ZetaRho,ZetaSigma,3,ZetaPts);
ZetaGrid = exp(ZetaGrid+ZetaMu);
%}
ZetaGrid = linspace(-3,3,ZetaPts);
% ZetaGrid = linspace(-2,2,5);
% ZetaPts = 5;
ZetaProb = tauchen_fg(0,ZetaGrid,0,1,1);
ZetaGrid = ZetaGrid*ZetaSigma + ZetaMu;
ZetaGrid = exp(ZetaGrid);
% attach entre exit shock
ZetaGrid = [0 ZetaGrid];
ZetaTrans = [
    (1-EntreShock1) EntreShock1*ZetaProb
    (1-EntreShock2)*ones(ZetaPts,1) EntreShock2*eye(ZetaPts)
    ];
ZetaPts = ZetaPts+1;
%}

EpsilonZetaTrans = kron(ZetaTrans,EpsilonTrans);

% wealth grid
AKsGrid = csvread('points.kgd');
AKsGrid = AKsGrid - AKsGrid(1) + 1e-1;

% AGrid = AKsGrid(1:6:end)/5;
AGrid = AKsGrid(1:10:end)/10;
APts = length(AGrid);
AMin = min(AGrid);
AMax = max(AGrid);

% tax policy
TaxBar = GBar + BBar*RBar;
Rho0 = TaxBar - Rho1*BBar - Rho2*GBar;

% wage
KLRatioBar = ((RBar+Delta)/Gamma)^(1/(Gamma-1));
WBar = (1-Gamma)*KLRatioBar^Gamma;

% VFI commonly used structure
[EpsilonMesh, ZetaMesh, AMesh] = ...
    ndgrid(EpsilonGrid, ZetaGrid, AGrid);
NX = numel(EpsilonMesh(:));
Worker.options.Algorithm = 'golden';
Worker.options.TolX = TolOpt;
Worker.options.ViolationLimit = 1e-6;
Worker.options.NumThreads = NumOfThreads;
Manager = Worker;
% asset goes last in memory
Idx = repmat([1:EpsilonPts*ZetaPts]',1,APts);
Idx = Idx(:)';
Worker.lb = repmat(AMin,1,NX);
Worker.idx = Idx;
Worker.rhs_return = ones(3,NX);
Manager.lb = [
    AMin*ones(1,NX);
    ];
Manager.idx = Idx;
Manager.rhs_return = ones(3,NX);

NewParams = v2struct;
end