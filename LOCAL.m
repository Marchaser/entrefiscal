function LOCAL(Params)
% solve vfi and simulation given steady state parameters
global SysEvalTimes;
SysEvalTimes = 0;
Params = COMMON(Params);
v2struct(Params);
global VfiRsltSs SmltRsltSs EVSs VWorkerSs VManagerSs zSs rSs wSs TaxSs LambdaSs DistSs GSs BSs

load('VfiRslt');
Params.TolOpt = eps;
Params.TolVfi = 1e-12;
Params.TolEqSs = eps;
Params.ShowDetail = 1;
VfiRsltSs = VFI_SS(ZBar,RBar,WBar,TaxBar,LambdaBar,TauLBar,TauRBar,TauPiBar,Params,VfiRslt.EV,[]);
SmltRsltSs = SIMULATE_SS(ZBar,BBar,GBar,VfiRsltSs,Params,[],[]);
Params.ShowDetail = 0;

EVSs = VfiRsltSs.EV;
VWorkerSs = VfiRsltSs.VWorker;
VManagerSs = VfiRsltSs.VManager;
zSs = ZBar;
DistSs = SmltRsltSs.Dist;
% DistSs(1) = 1-sum(DistSs(2:end));
rSs = SmltRsltSs.r;
wSs = SmltRsltSs.w;
TaxSs = TaxBar;
LambdaSs = LambdaBar;
GSs = GBar;
BSs = BBar;

Ks = SmltRsltSs.Ks;
Kd = SmltRsltSs.Kd;
K = SmltRsltSs.K;
Ns = SmltRsltSs.Ns;
Nd = SmltRsltSs.Nd;
N = SmltRsltSs.N;
r = SmltRsltSs.r;
w = SmltRsltSs.w;
Tax = TaxBar;
B = BBar;
z = ZBar;
G = GBar;
Lambda = LambdaBar;
IE = SmltRsltSs.IE;
IF = SmltRsltSs.IF;
I = SmltRsltSs.I;
YE = SmltRsltSs.YE;
YF = SmltRsltSs.YF;
Y = SmltRsltSs.Y;
EntrePopShare = SmltRsltSs.EntrePopShare;

Ss = [DistSs(:)' z G Lambda B];
Ss = [Ss EVSs(:)' VWorkerSs(:)' VManagerSs(:)'];
% Ss = [Ss Ks Kd K Ns Nd N r w];
Ss = [Ss Ks Kd K Ns Nd N r w Tax EntrePopShare];
nx = length(Ss);
nepsilon = 3;
neta = 2*EpsilonPts*ZetaPts*APts;

F = sys_stack([Ss Ss zeros(1,nepsilon) zeros(1,neta)],nx,nepsilon,neta,Params);
display(max(abs(F(:))));

Epsilon0 = zeros(1,nepsilon);
Eta0 = zeros(1,neta);

%{
% diagnostic what relative distance should be used
min(reshape(abs((VWorkerSs(:,:,2:end) - VWorkerSs(:,:,1:end-1))),1,[]))
min(reshape(abs((VManagerSs(:,:,2:end) - VManagerSs(:,:,1:end-1))),1,[]))
min(reshape(abs((EVSs(:,:,2:end) - EVSs(:,:,1:end-1))),1,[]))

min(reshape(abs((VWorkerSs(:,:,2:end) - VWorkerSs(:,:,1:end-1)) ./ VWorkerSs(:,:,1:end-1)),1,[]))
min(reshape(abs((VManagerSs(:,:,2:end) - VManagerSs(:,:,1:end-1)) ./ VManagerSs(:,:,1:end-1)),1,[]))
min(reshape(abs((EVSs(:,:,2:end) - EVSs(:,:,1:end-1)) ./ EVSs(:,:,1:end-1)),1,[]))
%}

% when using autodifferentiation, one should be careful preserving
% monotonicity
TypicalVWorker = VWorkerSs(:,:,2:end) - VWorkerSs(:,:,1:end-1);
TypicalVWorkerForward = cat(3,TypicalVWorker,inf*ones(EpsilonPts,ZetaPts,1));
TypicalVWorkerBackward = cat(3,inf*ones(EpsilonPts,ZetaPts,1),TypicalVWorker);
TypicalVWorker = min(TypicalVWorkerForward,TypicalVWorkerBackward);

TypicalVManager = VManagerSs(:,:,2:end) - VManagerSs(:,:,1:end-1);
TypicalVManagerForward = cat(3,TypicalVManager,inf*ones(EpsilonPts,ZetaPts,1));
TypicalVManagerBackward = cat(3,inf*ones(EpsilonPts,ZetaPts,1),TypicalVManager);
TypicalVManager = min(TypicalVManagerForward,TypicalVManagerBackward);

TypicalEV = EVSs(:,:,2:end) - EVSs(:,:,1:end-1);
TypicalEVForward = cat(3,TypicalEV,inf*ones(EpsilonPts,ZetaPts,1));
TypicalEVBackward = cat(3,inf*ones(EpsilonPts,ZetaPts,1),TypicalEV);
TypicalEV = min(TypicalEVForward,TypicalEVBackward);

TypicalX = [max(DistSs(:))*ones(1,EpsilonPts*ZetaPts*APts) z G Lambda B];
TypicalX = [TypicalX TypicalEV(:)' TypicalVWorker(:)' TypicalVManager(:)'];
% TypicalX = [TypicalX ones(1,3*EpsilonPts*ZetaPts*APts)];
TypicalX = [TypicalX Ks Kd K Ns Nd N r w Tax EntrePopShare];
TypicalX = [TypicalX TypicalX ones(1,nepsilon) ones(1,neta)];
% TypicalX = [ones(1,2*nx) ones(1,nepsilon) ones(1,neta)];
RelDeltaX = [1e-3*ones(1,EpsilonPts*ZetaPts*APts) 1e-3 1e-3 1e-3 1e-3];
RelDeltaX = [RelDeltaX 1e-3*ones(1,EpsilonPts*ZetaPts*APts) 1e-3*ones(1,EpsilonPts*ZetaPts*APts) 1e-3*ones(1,EpsilonPts*ZetaPts*APts)];
RelDeltaX = [RelDeltaX 1e-3*ones(1,10)];
RelDeltaX = [RelDeltaX RelDeltaX ones(1,nepsilon) ones(1,neta)];
assert(numel(TypicalX)==numel(RelDeltaX));
tic;
[F,Fp] = AutoDiff(@(x) sys_stack([x(1:nx) x(nx+1:2*nx) x(2*nx+1:2*nx+nepsilon) x(2*nx+nepsilon+1:2*nx+nepsilon+neta)], ...
    nx,nepsilon,neta,Params),[Ss Ss Epsilon0 Eta0],RelDeltaX,'central',TypicalX);
toc;

% call gensys
g0 = Fp(:,1:nx);
g1 = -Fp(:,nx+1:2*nx);
psi = -Fp(:,2*nx+1:2*nx+nepsilon);
pi = -Fp(:,2*nx+nepsilon+1:2*nx+nepsilon+neta);
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
IRF2(:,1) = IRF2(:,1) + impact(:,2);
for t=1:T
    IRF2(:,t+1) = G1*IRF2(:,t)+C;
end

IRF3 = IRF0;
IRF3(:,1) = IRF3(:,1) - impact(:,3);
for t=1:T
    IRF3(:,t+1) = G1*IRF3(:,t)+C;
end

% convert coefficients to canonical form
bp = EpsilonPts*ZetaPts*APts-1 + 4;
yidx = [bp+1:length(Ss)];
OmegaIdx = [1:bp];

D = C(yidx);
F = G1(yidx,OmegaIdx);

C1 = C(OmegaIdx);
A1 = G1(OmegaIdx,OmegaIdx);
A2 = G1(OmegaIdx,yidx);

C2 = C1+A2*D;
A = A1+A2*F;

SsOmega = (eye(size(A,1))-A)\C2;

save sys.mat Ss g0 g1 psi pi c F Fp G1 C impact fmat fwt ywt gev eu loose IRF1 IRF2 IRF3 D F A1 A2 C1 C2 A yidx OmegaIdx SsOmega
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

Dist = zeros(EpsilonPts,ZetaPts,APts);
Dist(:) = X(bp:bp+EpsilonPts*ZetaPts*APts-1);
% Dist(1) = 1-sum(Dist(2:end));
Distl = zeros(EpsilonPts,ZetaPts,APts);
Distl(:) = Xl(bp:bp+EpsilonPts*ZetaPts*APts-1);
% Distl(1) = 1-sum(Distl(2:end));
bp = bp + EpsilonPts*ZetaPts*APts;

z = X(bp);
zl = Xl(bp);
bp = bp+1;

G = X(bp);
Gl = Xl(bp);
bp = bp+1;

Lambda = X(bp);
Lambdal = Xl(bp);
bp = bp+1;

B = X(bp);
Bl = Xl(bp);
bp = bp+1;

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

Ks = X(bp);
Ksl = Xl(bp);
bp = bp+1;

Kd = X(bp);
Kdl = Xl(bp);
bp = bp+1;

K = X(bp);
Kl = Xl(bp);
bp = bp+1;

Ns = X(bp);
Nsl = Xl(bp);
bp = bp+1;

Nd = X(bp);
Ndl = Xl(bp);
bp = bp+1;

N = X(bp);
Nl = Xl(bp);
bp = bp+1;

r = X(bp);
rl = Xl(bp);
bp = bp+1;

w = X(bp);
wl = Xl(bp);
bp = bp+1;

Tax = X(bp);
Taxl = Xl(bp);
bp = bp+1;

EntrePopShare = X(bp);
EntrePopSharel = Xl(bp);
bp = bp+1;

assert(bp==length(X)+1);
% Taxl = Rho0 + Rho1*Bl + Rho2*Gl;

% VfiRslt = VFI_SS(zl,rl,wl,Taxl,Lambdal,TauLBar,TauRBar,TauPiBar,Params,EV,1);
global EVSs zSs rSs wSs TaxSs LambdaSs VfiRsltSs DistSs
if ~(isequal(EV,EVSs) && isequal(zl,zSs) && isequal(rl,rSs) && isequal(wl,wSs) && isequal(Taxl,TaxSs) && isequal(Lambdal,LambdaSs))
    % solve maximization problem
    VfiRslt = VFI_SS(zl,rl,wl,Taxl,Lambdal,TauLBar,TauRBar,TauPiBar,Params,EV,1);
else
    VfiRslt = VfiRsltSs;
end

%% update dist
% SmltRslt = SIMULATE_SS(zl,Bl,Gl,VfiRslt,Params,Distl,1);
global SmltRsltSs BSs GSs;
if ~(isequal(VfiRslt,VfiRsltSs) && isequal(Distl,DistSs) && isequal(zl,zSs) && isequal(Bl,BSs) && isequal(Gl,GSs))
    SmltRslt = SIMULATE_SS(zl,Bl,Gl,VfiRslt,Params,Distl,1);
else
    SmltRslt = SmltRsltSs;
end
%}

%% input to system of equations
bp = 1;

F(bp:bp+EpsilonPts*ZetaPts*APts-1) = reshape(SmltRslt.Dist-Dist,1,[]);
bp = bp+EpsilonPts*ZetaPts*APts;

F(bp) = log(z) - (ZRho*log(zl)+epsilon(1));
bp = bp+1;

F(bp) = G - (GRho*Gl+(1-GRho)*GBar+epsilon(2));
bp = bp+1;

F(bp) = Lambda - (LambdaRho*Lambdal+(1-LambdaRho)*LambdaBar+epsilon(3));
bp = bp+1;

F(bp) = B - (Bl*(1+rl)+Gl-Taxl);
bp = bp+1;

%{
T_EVl = VManagerl + normcdf((VWorkerl-VManagerl)/USigma).*(VWorkerl-VManagerl) ...
    +USigma*normpdf((VWorkerl-VManagerl)/USigma);
%}
T_EVl = VWorkerl + log(1 + exp(VManagerl-VWorkerl));
T_EVl = reshape(T_EVl,EpsilonPts,ZetaPts,APts);
T_EVl(:,1,:) = VWorkerl(:,1,:);
T_EVl = Beta*EpsilonZetaTrans * reshape(T_EVl, EpsilonPts*ZetaPts, APts);
T_EVl = reshape(T_EVl,EpsilonPts,ZetaPts,APts);
F(bp:bp+EpsilonPts*ZetaPts*APts-1) = reshape(T_EVl-EV,1,[]);
bp = bp+EpsilonPts*ZetaPts*APts;
%}

F(bp:bp+EpsilonPts*ZetaPts*APts-1) = reshape(VfiRslt.VWorker-VWorkerl,1,[])+eta(1:EpsilonPts*ZetaPts*APts);
bp = bp+EpsilonPts*ZetaPts*APts;

F(bp:bp+EpsilonPts*ZetaPts*APts-1) = reshape(VfiRslt.VManager-VManagerl,1,[])+eta(EpsilonPts*ZetaPts*APts+1:2*EpsilonPts*ZetaPts*APts);
bp = bp+EpsilonPts*ZetaPts*APts;

F(bp) = Ksl - SmltRslt.Ks;
bp = bp+1;

F(bp) = Kdl - SmltRslt.Kd;
bp = bp+1;

F(bp) = Kl - SmltRslt.K;
bp = bp+1;

F(bp) = Nsl - SmltRslt.Ns;
bp = bp+1;

F(bp) = Ndl - SmltRslt.Nd;
bp = bp+1;

F(bp) = Nl - SmltRslt.N;
bp = bp+1;

KLRatiol = Kl/Nl;
T_rl = zl*Gamma*KLRatiol^(Gamma-1)-Delta;
T_wl = zl*(1-Gamma)*KLRatiol^Gamma;
F(bp) = rl - T_rl;
bp = bp+1;

F(bp) = wl - T_wl;
bp = bp+1;

F(bp) = Taxl - (Rho0 + Rho1*Bl + Rho2*Gl);
bp = bp+1;

F(bp) = EntrePopSharel - SmltRslt.EntrePopShare;
bp = bp+1;
end
