function LOCAL(Params)
% solve vfi and simulation given steady state parameters
global SysEvalTimes;
SysEvalTimes = 0;
Params = COMMON(Params);
v2struct(Params);
global VfiRsltSs EVSs VWorkerSs VManagerSs zSs rSs wSs TaxSs LambdaSs DistSs

load('VfiRslt');
VfiRsltSs = VFI_SS(ZBar,RBar,WBar,TaxBar,LambdaBar,TauLBar,TauRBar,TauPiBar,Params,VfiRslt.EV,[]);
SmltRslt = SIMULATE_SS(BBar,GBar,VfiRsltSs,Params,[],[]);

EVSs = VfiRsltSs.EV;
VWorkerSs = VfiRsltSs.VWorker;
VManagerSs = VfiRsltSs.VManager;
zSs = ZBar;
DistSs = SmltRslt.Dist;
rSs = SmltRslt.r;
wSs = SmltRslt.w;
TaxSs = TaxBar;
LambdaSs = LambdaBar;

Ks = SmltRslt.Ks;
Kd = SmltRslt.Kd;
K = SmltRslt.K;
Ns = SmltRslt.Ns;
Nd = SmltRslt.Nd;
N = SmltRslt.N;
r = SmltRslt.r;
w = SmltRslt.w;
Tax = TaxBar;
B = BBar;
z = ZBar;
G = GBar;
Lambda = LambdaBar;
IE = SmltRslt.IE;
IF = SmltRslt.IF;
I = SmltRslt.I;
YE = SmltRslt.YE;
YF = SmltRslt.YF;
Y = SmltRslt.Y;
EntrePopShare = SmltRslt.EntrePopShare;

SsOthers = [Ks Kd K Ns Nd N r w Tax B z G Lambda];
SsOthers = [SsOthers YE YF Y EntrePopShare];
% SsOthers = log(SsOthers);
Ss = [EVSs(:)' VWorkerSs(:)' VManagerSs(:)' DistSs(:)' SsOthers];
nx = length(Ss);
nepsilon = 3;
neta = 2*EpsilonPts*ZetaPts*APts;

F = sys_stack([Ss Ss zeros(1,nepsilon) zeros(1,neta)],nx,nepsilon,neta,Params);
display(max(abs(F(:))));

Epsilon0 = zeros(1,nepsilon);
Eta0 = zeros(1,neta);
tic;

TypicalX = max(abs([Ss Ss ones(1,nepsilon) ones(1,neta)]),1e-2);
[F,Fp] = AutoDiff(@(x) sys_stack([x(1:nx) x(nx+1:2*nx) x(2*nx+1:2*nx+nepsilon) x(2*nx+nepsilon+1:2*nx+nepsilon+neta)], ...
    nx,nepsilon,neta,Params),[Ss Ss Epsilon0 Eta0],1e-5,'central',TypicalX);
toc;

% call gensys
g0 = Fp(:,1:nx);
g1 = -Fp(:,nx+1:2*nx);
psi = -Fp(:,2*nx+1:2*nx+nepsilon);
pi = -Fp(:,2*nx+nepsilon+1:end);
c = (Fp(:,1:nx)+Fp(:,nx+1:2*nx))*Ss';

t2 = tic;
[G1,C,impact,fmat,fwt,ywt,gev,eu,loose] = gensys(g0,g1,c,psi,pi);
toc(t2);

T = 100;
IRF0 = repmat(Ss(:),1,T+1);

IRF1 = IRF0;
IRF1(:,1) = IRF1(:,1) + impact(:,1);
for t=1:T
    IRF1(:,t+1) = G1*IRF1(:,t)+C;
end

IRF2 = IRF0;
IRF2(:,9) = IRF2(:,9) + impact(:,2);
for t=9:T
    IRF2(:,t+1) = G1*IRF2(:,t)+C;
end

IRF3 = IRF0;
IRF3(:,1) = IRF3(:,1) - impact(:,3);
for t=1:T
    IRF3(:,t+1) = G1*IRF3(:,t)+C;
end

IRFLambda1 = IRF3;
IRFLambda1(:,9) = IRF3(:,9) + impact(:,2);
for t=9:T
    IRFLambda1(:,t+1) = G1*IRFLambda1(:,t)+C;
end

save sys.mat Ss g0 g1 psi pi c F Fp G1 C impact fmat fwt ywt gev eu loose IRF1 IRF2 IRF3
end

function FStack = sys_stack(x,nx,nepsilon,neta,Params)
X = x(1:nx);
Xl = x(nx+1:2*nx);
epsilon = x(2*nx+1:2*nx+nepsilon);
eta = x(2*nx+nepsilon+1:2*nx+nepsilon+neta);
FStack = sys(X,Xl,epsilon,eta,Params);
end

function F = sys(X,Xl,epsilon,eta,Params)
global SysEvalTimes
SysEvalTimes = SysEvalTimes+1;
if mod(SysEvalTimes,100) == 1
    display(SysEvalTimes);
end
v2struct(Params);

bp = 1;
EV = X(bp:bp+EpsilonPts*ZetaPts*APts-1);
EV = reshape(EV,EpsilonPts,ZetaPts,APts);
EVl = Xl(bp:bp+EpsilonPts*ZetaPts*APts-1);
EVl = reshape(EVl,EpsilonPts,ZetaPts,APts);
bp = bp + EpsilonPts*ZetaPts*APts;

VWorker = X(bp:bp+EpsilonPts*ZetaPts*APts-1);
VWorker = reshape(VWorker,EpsilonPts,ZetaPts,APts);
VWorkerl = Xl(bp:bp+EpsilonPts*ZetaPts*APts-1);
VWorkerl = reshape(VWorkerl,EpsilonPts,ZetaPts,APts);
bp = bp + EpsilonPts*ZetaPts*APts;

VManager = X(bp:bp+EpsilonPts*ZetaPts*APts-1);
VManager = reshape(VManager,EpsilonPts,ZetaPts,APts);
VManagerl = Xl(bp:bp+EpsilonPts*ZetaPts*APts-1);
VManagerl = reshape(VManagerl,EpsilonPts,ZetaPts,APts);
bp = bp + EpsilonPts*ZetaPts*APts;

Dist = X(bp:bp+EpsilonPts*ZetaPts*APts-1);
Dist = reshape(Dist,[EpsilonPts ZetaPts APts]);
Distl = Xl(bp:bp+EpsilonPts*ZetaPts*APts-1);
Distl = reshape(Distl,[EpsilonPts ZetaPts APts]);
bp = bp + EpsilonPts*ZetaPts*APts;

%{
% log from here below
X(bp:end) = exp(X(bp:end));
Xl(bp:end) = exp(Xl(bp:end));
%}

Ks = X(bp);
Ksl = Xl(bp);
bp = bp+1;

%19
Kd = X(bp);
Kdl = Xl(bp);
bp = bp+1;

%18
K = X(bp);
Kl = Xl(bp);
bp = bp+1;

%17
Ns = X(bp);
Nsl = Xl(bp);
bp = bp+1;

%16
Nd = X(bp);
Ndl = Xl(bp);
bp = bp+1;

%15
N = X(bp);
Nl = Xl(bp);
bp = bp+1;

%14
r = X(bp);
rl = Xl(bp);
bp = bp+1;

%13
w = X(bp);
wl = Xl(bp);
bp = bp+1;

%12
Tax = X(bp);
Taxl = Xl(bp);
bp = bp+1;

%11
B = X(bp);
Bl = Xl(bp);
bp = bp+1;

%10
z = X(bp);
zl = Xl(bp);
bp = bp+1;

%9
G = X(bp);
Gl = Xl(bp);
bp = bp+1;

%8
Lambda = X(bp);
Lambdal = Xl(bp);
bp = bp+1;

%4
YE = X(bp);
YEl = Xl(bp);
bp = bp+1;

%3
YF = X(bp);
YFl = Xl(bp);
bp = bp+1;

%2
Y = X(bp);
Yl = Xl(bp);
bp = bp+1;

%1
EntrePopShare = X(bp);
EntrePopSharel = Xl(bp);
bp = bp+1;

assert(bp==length(X)+1);

% EV = Beta*EpsilonZetaTrans * reshape(V, EpsilonPts*ZetaPts, APts);
% VfiResult = VFI_SS(zl,rl,wl,Taxl,Lambdal,TauLBar,TauRBar,TauPiBar,Params,EV,1);
global EVSs zSs rSs wSs TaxSs LambdaSs VfiRsltSs DistSs
if ~(isequal(EV,EVSs) && isequal(zl,zSs) && isequal(rl,rSs) && isequal(wl,wSs) && isequal(Taxl,TaxSs) && isequal(Lambdal,LambdaSs))
    %% solve maximization problem
    % EV = Beta*EpsilonZetaTrans * reshape(V, EpsilonPts*ZetaPts, APts);
    VfiRslt = VFI_SS(zl,rl,wl,Taxl,Lambdal,TauLBar,TauRBar,TauPiBar,Params,EV,1);
else
    VfiRslt = VfiRsltSs;
end
%}

T_VWorkerl = VfiRslt.VWorker;
T_VManagerl = VfiRslt.VManager;
OccPolicy = VfiRslt.OccPolicy;

%% update dist
FullPp = tensor_pchip({AGrid}, cat(1,...
    reshape(VfiRslt.ApWorker, 1, EpsilonPts, ZetaPts, APts), ...
    reshape(VfiRslt.ApManager, 1, EpsilonPts, ZetaPts, APts)));
FullPp = myppual(FullPp);

ApInterp = myppualMKL_CMEX(int32(NumOfThreads), {AGrid}, FullPp.coefs, [], int32([4]), int32(2*EpsilonPts*ZetaPts), [], [AGrid(:)'], [], [], []);
ApInterp = reshape(ApInterp, 2, EpsilonPts, ZetaPts, APts);
ApInterp = min(max(ApInterp,AMin),AMax);
[~,ApCell] = histc(ApInterp,AGrid);
ApCell = min(ApCell,APts-1);
ApLeftShare = (AGrid(ApCell+1)-ApInterp) ./ (AGrid(ApCell+1)-AGrid(ApCell));

ApWorkerCell = reshape(ApCell(1,:),EpsilonPts,ZetaPts,APts);
ApManagerCell = reshape(ApCell(2,:),EpsilonPts,ZetaPts,APts);
ApWorkerLeftShare = reshape(ApLeftShare(1,:),EpsilonPts,ZetaPts,APts);
ApManagerLeftShare = reshape(ApLeftShare(2,:),EpsilonPts,ZetaPts,APts);

ApWorkerCellInt = int32(ApWorkerCell)-1;
ApManagerCellInt = int32(ApManagerCell)-1;

% T_Dist = update_dist(Distl,OccPolicy,ApManagerCellInt,ApWorkerCellInt,ApManagerLeftShare,ApWorkerLeftShare,EpsilonTrans,ZetaTrans);
if ~(isequal(VfiRslt,VfiRsltSs) && isequal(Distl,DistSs))
    T_Dist = update_dist(Distl,OccPolicy,ApManagerCellInt,ApWorkerCellInt,ApManagerLeftShare,ApWorkerLeftShare,EpsilonTrans,ZetaTrans);
else
    T_Dist = DistSs;
end
%}

%% input to system of equations
bp = 1;
T_EVl = VManagerl + normcdf((VWorkerl-VManagerl)/USigma).*(VWorkerl-VManagerl) ...
    +USigma*normpdf((VWorkerl-VManagerl)/USigma);
T_EVl = reshape(T_EVl,EpsilonPts,ZetaPts,APts);
T_EVl(:,1,:) = VWorkerl(:,1,:);
T_EVl = Beta*EpsilonZetaTrans * reshape(T_EVl, EpsilonPts*ZetaPts, APts);
T_EVl = reshape(T_EVl,EpsilonPts,ZetaPts,APts);
F(bp:bp+EpsilonPts*ZetaPts*APts-1) = reshape(T_EVl-EVl,1,[]);
bp = bp+EpsilonPts*ZetaPts*APts;

F(bp:bp+EpsilonPts*ZetaPts*APts-1) = reshape(T_VWorkerl-VWorkerl,1,[])+eta(1:EpsilonPts*ZetaPts*APts);
bp = bp+EpsilonPts*ZetaPts*APts;

F(bp:bp+EpsilonPts*ZetaPts*APts-1) = reshape(T_VManagerl-VManagerl,1,[])+eta(EpsilonPts*ZetaPts*APts+1:2*EpsilonPts*ZetaPts*APts);
bp = bp+EpsilonPts*ZetaPts*APts;

F(bp:bp+EpsilonPts*ZetaPts*APts-1) = reshape(T_Dist-Dist,1,[]);
bp = bp+EpsilonPts*ZetaPts*APts;

%19
F(bp) = Ksl - sum(Distl(:).*AMesh(:));
bp = bp+1;

%18
F(bp) = Kdl - sum(Distl(:).*OccPolicy(:).*VfiRslt.KManager(:));
bp = bp+1;

%17
F(bp) = Kl - (Ksl-Bl-Kdl);
bp = bp+1;

%16
F(bp) = Nsl - sum(Distl(:).*(1-OccPolicy(:)).*VfiRslt.NWorker(:).*EpsilonMesh(:));
bp = bp+1;

%15
F(bp) = Ndl - sum(Distl(:).*OccPolicy(:).*VfiRslt.NManager(:));
bp = bp+1;

%14
F(bp) = Nl - (Nsl-Ndl);
bp = bp+1;

%13
KLRatiol = Kl/Nl;
F(bp) = rl - (zl*Gamma*KLRatiol^(Gamma-1)-Delta);
bp = bp+1;

%12
F(bp) = wl - zl*(1-Gamma)*KLRatiol^Gamma;
bp = bp+1;

%11
F(bp) = Taxl - (Rho0+Rho1*Bl+Rho2*Gl);
bp = bp+1;

%10
F(bp) = B - (Bl*(1+rl)+Gl-Taxl);
bp = bp+1;

%9
F(bp) = log(z) - (ZRho*log(zl)+epsilon(1));
bp = bp+1;

%8
F(bp) = G - (GRho*Gl+(1-GRho)*GBar+epsilon(2));
bp = bp+1;

%7
F(bp) = Lambda - (LambdaRho*Lambdal+(1-LambdaRho)*LambdaBar+epsilon(3));
bp = bp+1;

%3
F(bp) = YEl - sum(Distl(:).*VfiRslt.YManager(:).*OccPolicy(:));
bp = bp+1;

%2
F(bp) = YFl - zl*Kl^Gamma*Nl^(1-Gamma);
bp = bp+1;

%1
F(bp) = Yl - (YEl+YFl);
bp = bp+1;

%
F(bp) = EntrePopSharel - sum(Distl(:).*VfiRslt.OccPolicy(:));
bp = bp+1;
end
