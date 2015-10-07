function SIMULATE(sys,T,shock,Params)
Params = COMMON(Params);
% Params.ShowDetail = 1;

Ss = sys.Ss';
G1 = sys.G1;
C = sys.C;
%{
for j=1:100
    Ss = G1*Ss + C;
end
%}
impact = sys.impact;

Beta = Params.Beta;
EpsilonPts = Params.EpsilonPts;
ZetaPts = Params.ZetaPts;
APts = Params.APts;
EpsilonZetaTrans = Params.EpsilonZetaTrans;
Gamma = Params.Gamma;
Delta = Params.Delta;
TauLBar = Params.TauLBar;
TauRBar = Params.TauRBar;
TauPiBar = Params.TauPiBar;
Rho0 = Params.Rho0;
Rho1 = Params.Rho1;
Rho2 = Params.Rho2;
KLRatioMin = ((1/Beta-1 + Delta + 1e-3)/Gamma)^(1/(Gamma-1));
KLRatioMax = (Delta/Gamma)^(1/(Gamma-1));
SmltRslt = [];
ZRho = Params.ZRho;
GRho = Params.GRho;
LambdaRho = Params.LambdaRho;
ZBar = Params.ZBar;
GBar = Params.GBar;
LambdaBar = Params.LambdaBar;
USigma = Params.USigma;

    function Metric = getKLRatioMetric(KLRatio)
        R = Gamma*KLRatio^(Gamma-1)-Delta;
        W = (1-Gamma)*KLRatio^Gamma;
        % solve decision problem
        VfiRslt = VFI_SS(z(t),R,W,Tax(t),Lambda(t),TauLBar,TauRBar,TauPiBar,Params,EV,1);
        % simulate
        SmltRslt = SIMULATE_SS(z(t),B(t),G(t),VfiRslt,Params,Dist,1);
        % compute metric
        Metric = SmltRslt.KLRatio - KLRatio;
    end

Dist = zeros(EpsilonPts,ZetaPts,APts);
% Dist(:) = Ss(1:EpsilonPts*ZetaPts*APts);
% Ss(1) = 1 - sum(Dist(2:end));
Response(:,1) = Ss;
for t=1:T
    % get states
    bp = 1;
    
    Dist(:) = Response(bp:bp-1+EpsilonPts*ZetaPts*APts,t);
    bp = bp + EpsilonPts*ZetaPts*APts;
    
    z(t) = Response(bp,t);
    bp = bp+1;
    
    G(t) = Response(bp,t);
    bp = bp+1;
    
    Lambda(t) = Response(bp,t);
    bp = bp+1;
    
    B(t) = Response(bp,t);
    bp = bp+1;
    
    % other endogenous varialbe directly computed
    Tax(t) = Rho0+Rho1*B(t)+Rho2*G(t);

    bbp = bp; % base pointer from mother
    % predict future values
    Response(:,t+1) = G1*Response(:,t) + C;
    
    VWorker = Response(bbp:bbp-1+EpsilonPts*ZetaPts*APts,t+1);
    VWorker = reshape(VWorker,EpsilonPts,ZetaPts,APts);
    bbp = bbp+EpsilonPts*ZetaPts*APts;
    
    VManager = Response(bbp:bbp-1+EpsilonPts*ZetaPts*APts,t+1);
    VManager = reshape(VManager,EpsilonPts,ZetaPts,APts);
    bbp = bbp+EpsilonPts*ZetaPts*APts;
    
    r(t) = Response(bbp,t);
    bbp = bbp+1;
    w(t) = Response(bbp,t);
    bbp = bbp+1;
    
    EV = VManager + normcdf((VWorker-VManager)/USigma).*(VWorker-VManager) ...
        +USigma*normpdf((VWorker-VManager)/USigma);
    EV = reshape(EV,EpsilonPts,ZetaPts,APts);
    EV(:,1,:) = VWorker(:,1,:);
    EV = Beta*EpsilonZetaTrans * reshape(EV, EpsilonPts*ZetaPts, APts);
    EV = reshape(EV,EpsilonPts,ZetaPts,APts);
    
    % find KLRatio that clears market
    % options = optimoptions('lsqnonlin','Display','iter','TolX',1e-8,'TolFun',1e-20);
    options = optimset('Display','iter','TolX',eps);
    KLRatio = fzero(@getKLRatioMetric,[KLRatioMin,KLRatioMax],options);
    %}
    %{
    VfiRslt = VFI_SS(z(t),r(t),w(t),Tax(t),Lambda(t),TauLBar,TauRBar,TauPiBar,Params,EV,1);
    SmltRslt = SIMULATE_SS(z(t),B(t),G(t),VfiRslt,Params,Dist,1);
    %}
    
    % jump varialbe shoudln't mater. If it does, specify here.
    % get aggregates
    Ks(t) = SmltRslt.Ks;
    Kd(t) = SmltRslt.Kd;
    K(t) = SmltRslt.K;
    Ns(t) = SmltRslt.Ns;
    Nd(t) = SmltRslt.Nd;
    N(t) = SmltRslt.N;
    r(t) = SmltRslt.r;
    w(t) = SmltRslt.w;
    YE(t) = SmltRslt.YE;
    YF(t) = SmltRslt.YF;
    Y(t) = SmltRslt.Y;
    EntrePopShare(t) = SmltRslt.EntrePopShare;
    
    % updates current jump variables
    bbp = bp;
    Response(bbp:bbp-1+EpsilonPts*ZetaPts*APts,t) = VfiRslt.VWorker(:);
    bbp = bbp+EpsilonPts*ZetaPts*APts;
    
    Response(bbp:bbp-1+EpsilonPts*ZetaPts*APts,t) = VfiRslt.VManager(:);
    bbp = bbp+EpsilonPts*ZetaPts*APts;
    
    Response(bbp,t) = SmltRslt.r;
    bbp = bbp+1;
    
    Response(bbp,t) = SmltRslt.w;
    bbp = bbp+1;
    %}
    % evolution of states, inexact
    Response(:,t+1) = G1*Response(:,t) + C + impact*shock(:,t);
    
    % evolution of states, exact
    bp = 1;
    Response(bp:bp-1+EpsilonPts*ZetaPts*APts,t+1) = SmltRslt.Dist(:);
    bp = bp+EpsilonPts*ZetaPts*APts;
    
    Response(bp,t+1) = exp(log(z(t))*ZRho + (1-ZRho)*log(ZBar) + shock(1,t));
    bp = bp+1;
    
    Response(bp,t+1) = GRho*G(t) + (1-GRho)*GBar + shock(2,t);
    bp = bp+1;
    
    Response(bp,t+1) = LambdaRho*Lambda(t) + (1-LambdaRho)*LambdaBar + shock(3,t);
    bp = bp+1;
    
    Response(bp,t+1) = SmltRslt.Bp;
    bp = bp+1;
    %}
end
end