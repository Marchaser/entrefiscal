function LOCAL(Params)
% solve vfi and simulation given steady state parameters
global SysEvalTimes;
SysEvalTimes = 0;
Params = COMMON(Params);
v2struct(Params);
global VfiRsltSs VSs zSs rSs wSs TaxSs LambdaSs DistSs

load('VfiRslt');
VfiRsltSs = VFI_SS(ZBar,RBar,WBar,TaxBar,LambdaBar,TauLBar,TauRBar,TauPiBar,Params,VfiRslt.EV,[]);
SmltRslt = SIMULATE_SS(BBar,GBar,VfiRsltSs,Params,[],[]);

VSs = VfiRsltSs.V;
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
SsOthers = [SsOthers IE IF I];
SsOthers = [SsOthers YE YF Y EntrePopShare];
SsOthers = log(SsOthers);
Ss = [VSs(:)' DistSs(:)' SsOthers];
nx = length(Ss);
nepsilon = 3;
neta = EpsilonPts*ZetaPts*APts;

F = sys_stack([Ss Ss zeros(1,nepsilon) zeros(1,neta)],nx,nepsilon,neta,Params);
display(max(abs(F(:))));

Epsilon0 = zeros(1,nepsilon);
Eta0 = zeros(1,neta);
tic;
TypicalXSsOthers = abs(SsOthers);
TypicalXSsOthers(11) = 1e-2;
TypicalX = [VSs(:)' max(DistSs(:)')*ones(1,length(DistSs(:)')) TypicalXSsOthers ...
    VSs(:)' max(DistSs(:)')*ones(1,length(DistSs(:)')) TypicalXSsOthers ...
    ZBar GBar LambdaBar VSs(:)'];
RelDeltaX = [1e-5*ones(1,length(VSs(:)')) 1e-3*ones(1,length(DistSs(:)')) 1e-3*ones(1,length(SsOthers)) ...
    1e-5*ones(1,length(VSs(:)')) 1e-3*ones(1,length(DistSs(:)')) 1e-3*ones(1,length(SsOthers)) ...
    1 1 1 ones(1,length(VSs(:)'))];
[F,Fp] = AutoDiff(@(x) sys_stack([x(1:nx) x(nx+1:2*nx) x(2*nx+1:2*nx+nepsilon) x(2*nx+nepsilon+1:2*nx+nepsilon+neta)], ...
    nx,nepsilon,neta,Params),[Ss Ss Epsilon0 Eta0],RelDeltaX,'central',TypicalX);
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

IRF1(:,1) = impact(:,1);
for t=1:100
    IRF1(:,t+1) = G1*IRF1(:,t);
end

IRF2(:,1) = impact(:,2);
for t=1:100
    IRF2(:,t+1) = G1*IRF2(:,t);
end

IRF3(:,1) = impact(:,3);
for t=1:100
    IRF3(:,t+1) = G1*IRF3(:,t);
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
V = X(bp:bp+EpsilonPts*ZetaPts*APts-1);
V = reshape(V,EpsilonPts,ZetaPts,APts);
Vl = Xl(bp:bp+EpsilonPts*ZetaPts*APts-1);
bp = bp + EpsilonPts*ZetaPts*APts;
Vl = reshape(Vl,EpsilonPts,ZetaPts,APts);

Dist = X(bp:bp+EpsilonPts*ZetaPts*APts-1);
Dist = reshape(Dist,[EpsilonPts ZetaPts APts]);
Distl = Xl(bp:bp+EpsilonPts*ZetaPts*APts-1);
Distl = reshape(Distl,[EpsilonPts ZetaPts APts]);
bp = bp + EpsilonPts*ZetaPts*APts;

Ks = exp(X(bp));
Ksl = exp(Xl(bp));
bp = bp+1;

Kd = exp(X(bp));
Kdl = exp(Xl(bp));
bp = bp+1;

K = exp(X(bp));
Kl = exp(Xl(bp));
bp = bp+1;

Ns = exp(X(bp));
Nsl = exp(Xl(bp));
bp = bp+1;

Nd = exp(X(bp));
Ndl = exp(Xl(bp));
bp = bp+1;

N = exp(X(bp));
Nl = exp(Xl(bp));
bp = bp+1;

r = exp(X(bp));
rl = exp(Xl(bp));
bp = bp+1;

w = exp(X(bp));
wl = exp(Xl(bp));
bp = bp+1;

Tax = exp(X(bp));
Taxl = exp(Xl(bp));
bp = bp+1;

B = exp(X(bp));
Bl = exp(Xl(bp));
bp = bp+1;

z = exp(X(bp));
zl = exp(Xl(bp));
bp = bp+1;

G = exp(X(bp));
Gl = exp(Xl(bp));
bp = bp+1;

Lambda = exp(X(bp));
Lambdal = exp(Xl(bp));
bp = bp+1;

IE = exp(X(bp));
IEl = exp(Xl(bp));
bp = bp+1;

IF = exp(X(bp));
IFl = exp(Xl(bp));
bp = bp+1;

I = exp(X(bp));
Il = exp(Xl(bp));
bp = bp+1;

YE = exp(X(bp));
YEl = exp(Xl(bp));
bp = bp+1;

YF = exp(X(bp));
YFl = exp(Xl(bp));
bp = bp+1;

Y = exp(X(bp));
Yl = exp(Xl(bp));
bp = bp+1;

EntrePopShare = exp(X(bp));
EntrePopSharel = exp(Xl(bp));
bp = bp+1;

assert(bp==length(X)+1);

% EV = Beta*EpsilonZetaTrans * reshape(V, EpsilonPts*ZetaPts, APts);
% VfiResult = VFI_SS(zl,rl,wl,Taxl,Lambdal,TauLBar,TauRBar,TauPiBar,Params,EV,1);
global VSs zSs rSs wSs TaxSs LambdaSs VfiRsltSs DistSs
if ~(isequal(V,VSs) && isequal(zl,zSs) && isequal(rl,rSs) && isequal(wl,wSs) && isequal(Taxl,TaxSs) && isequal(Lambdal,LambdaSs))
    %% solve maximization problem
    EV = Beta*EpsilonZetaTrans * reshape(V, EpsilonPts*ZetaPts, APts);
    VfiResult = VFI_SS(zl,rl,wl,Taxl,Lambdal,TauLBar,TauRBar,TauPiBar,Params,EV,1);
else
    VfiResult = VfiRsltSs;
end
%}

T_Vl = VfiResult.V;
OccPolicy = VfiResult.OccPolicy;

%% update dist
FullPp = tensor_pchip({AGrid}, cat(1,...
    reshape(VfiResult.ApWorker, 1, EpsilonPts, ZetaPts, APts), ...
    reshape(VfiResult.ApManager, 1, EpsilonPts, ZetaPts, APts)));
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
if ~(isequal(VfiResult,VfiRsltSs) && isequal(Distl,DistSs)) && 0
    T_Dist = update_dist(Distl,OccPolicy,ApManagerCellInt,ApWorkerCellInt,ApManagerLeftShare,ApWorkerLeftShare,EpsilonTrans,ZetaTrans);
else
    T_Dist = DistSs;
end
%}

%% input to system of equations
bp = 1;
F(bp:bp+EpsilonPts*ZetaPts*APts-1) = reshape(T_Vl-Vl,1,[]) + eta;
bp = bp+EpsilonPts*ZetaPts*APts;

F(bp:bp+EpsilonPts*ZetaPts*APts-1) = reshape(T_Dist-Dist,1,[]);
bp = bp+EpsilonPts*ZetaPts*APts;

%19
F(bp) = Ksl - sum(Distl(:).*AMesh(:));
bp = bp+1;

%18
F(bp) = Kdl - sum(Distl(:).*OccPolicy(:).*VfiResult.KManager(:));
bp = bp+1;

%17
F(bp) = Kl - (Ksl-Bl-Kdl);
bp = bp+1;

%16
F(bp) = Nsl - sum(Distl(:).*(1-OccPolicy(:)).*VfiResult.NWorker(:).*EpsilonMesh(:));
bp = bp+1;

%15
F(bp) = Ndl - sum(Distl(:).*OccPolicy(:).*VfiResult.NManager(:));
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

%6
F(bp) = IE - (Kd-(1-Delta)*Kdl);
bp = bp+1;

%5
F(bp) = IF - (K-(1-Delta)*Kl);
bp = bp+1;

%4
F(bp) = I - (IE+IF);
bp = bp+1;

%3
F(bp) = YEl - sum(Distl(:).*VfiResult.YManager(:).*OccPolicy(:));
bp = bp+1;

%2
F(bp) = YFl - zl*Kl^Gamma*Nl^(1-Gamma);
bp = bp+1;

%1
F(bp) = Yl - (YEl+YFl);
bp = bp+1;

%
F(bp) = EntrePopSharel - sum(Distl(:).*VfiResult.OccPolicy(:));
bp = bp+1;
end