function RA
% parameters
Params = SETUP;
Params = COMMON(Params);
global Beta Sigma Chi Gamma Delta ZRho GRho Rho0 Rho1 Rho2 GBar
Beta = Params.Beta;
Sigma = Params.Sigma;
Chi = Params.Chi;
Gamma = Params.Gamma;
Delta = Params.Delta;
ZRho = Params.ZRho;
GRho = Params.GRho;
Rho0 = Params.Rho0;
Rho1 = Params.Rho1;
Rho2 = Params.Rho2;
GBar = Params.GBar;

Rho0 = 0.144;

% solve steady state
X0 = ones(1,12);
options = optimset('MaxFunEvals',inf,'MaxIter',inf,'Display','iter','FinDiffRelStep',1e-4,'TolFun',1e-20,'TolX',1e-20);

% [X,F] = lsqnonlin(@ss,X0,[0 0 0],[],options);

% compute Fp
nx = length(X0);
nepsilon = 2;
neta = 1;
[Ss,F] = fsolve(@(x) ss(x,nepsilon,neta),X0,options);

epsilon0 = zeros(1,nepsilon);
eta0 = zeros(1,neta);
[~,Fp] = AutoDiff(@(x) sys_stack([x(1:nx) x(nx+1:2*nx) x(2*nx+1:2*nx+nepsilon) x(2*nx+nepsilon+1:2*nx+nepsilon+neta)],nx,nepsilon,neta), [Ss Ss epsilon0 eta0], ...
    1e-6, 'central',ones(1,2*nx+nepsilon+neta));

% call gensys
g0 = Fp(:,1:nx);
g1 = -Fp(:,nx+1:2*nx);
psi = -Fp(:,2*nx+1:2*nx+nepsilon);
pi = -Fp(:,2*nx+nepsilon+1:2*nx+nepsilon+neta);
c = (Fp(:,1:nx)+Fp(:,nx+1:2*nx))*Ss';

[G1,C,impact,fmat,fwt,ywt,gev,eu] = gensys(g0,g1,c,psi,pi);



% compute impulse response
IRF(:,2) = impact(:,2);
for t=2:30
   IRF(:,t+1) = G1*IRF(:,t);
end

plot(IRF(end,:));
%}

save sysRA.mat G1 C impact fmat fwt ywt gev eu IRF
end

function F = ss(X,nepsilon,neta)
F = sys(X,X,zeros(1,nepsilon),zeros(1,neta));
end

function F = sys(X,Xl,Epsilon,Eta)
global Beta Sigma Chi Gamma Delta ZRho GRho Rho0 Rho1 Rho2 GBar

bp = 1;

z = X(bp);
zl = Xl(bp);
bp = bp+1;

G = X(bp);
Gl = Xl(bp);
bp = bp+1;

a = X(bp);
al = Xl(bp);
bp = bp+1;

B = X(bp);
Bl = Xl(bp);
bp = bp+1;

w = X(bp);
wl = Xl(bp);
bp = bp+1;

r = X(bp);
rl = Xl(bp);
bp = bp+1;

c = X(bp);
cl = Xl(bp);
bp = bp+1;

K = X(bp);
Kl = Xl(bp);
bp = bp+1;

N = X(bp);
Nl = Xl(bp);
bp = bp+1;

E = X(bp);
El = Xl(bp);
bp = bp+1;

Tax = X(bp);
Taxl = Xl(bp);
bp = bp+1;

Y = X(bp);
Yl = Xl(bp);
bp = bp+1;

assert(bp==length(X)+1);

bp=1;

F(bp) = log(z) - (ZRho*log(zl)+Epsilon(1));
bp = bp+1;

F(bp) = G - (GRho*Gl+(1-GRho)*GBar+Epsilon(2));
bp = bp+1;

F(bp) = B - (Bl*(1+rl)+Gl-Taxl);
bp = bp+1;

F(bp) = El - (wl+(1+rl)*al-Taxl-a);
bp = bp+1;

F(bp) = cl - Chi*El;
bp = bp+1;

F(bp) = Nl - (1-(1-Chi)*El/wl);
bp = bp+1;

F(bp) = al - (Kl+Bl);
bp = bp+1;

F(bp) = wl - (zl*(1-Gamma)*(Kl/Nl)^Gamma);
bp = bp+1;

F(bp) = rl - (zl*Gamma*(Kl/Nl)^(Gamma-1)-Delta);
bp = bp+1;

F(bp) = Yl - (zl*Kl^Gamma*Nl^(1-Gamma));
bp = bp+1;

F(bp) = Taxl - (Rho0+Rho1*Bl+Rho2*Gl);
bp = bp+1;

F(bp) = El^(-Sigma) - Beta*(1+r)*E^(-Sigma) + Eta;
bp = bp+1;

assert(bp==length(X)+1);
end

function FStack = sys_stack(x,nx,nepsilon,neta)
X = x(1:nx);
Xl = x(nx+1:2*nx);
epsilon = x(2*nx+1:2*nx+nepsilon);
eta = x(2*nx+nepsilon+1:2*nx+nepsilon+neta);
FStack = sys(X,Xl,epsilon,eta);
end

