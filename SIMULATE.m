function SimulateResult = SIMULATE(sys,T,shock,Params)
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
KLRatioLast = Params.KLRatioBar;

    function Metric = getKLRatioMetric(KLRatio)
        R = z(t)*Gamma*KLRatio^(Gamma-1)-Delta;
        W = z(t)*(1-Gamma)*KLRatio^Gamma;
        % solve decision problem
        VfiRslt = VFI_SS(z(t),R,W,Tax(t),Lambda(t),TauLBar,TauRBar,TauPiBar,Params,EV,1);
        % simulate
        SmltRslt = SIMULATE_SS(z(t),B(t),G(t),VfiRslt,Params,Dist,1);
        % compute metric
        Metric = SmltRslt.KLRatio - KLRatio;
    end

Dist = zeros(EpsilonPts,ZetaPts,APts);
Dist(:) = Ss(1:EpsilonPts*ZetaPts*APts);
Ss(1) = 1 - sum(Dist(2:end));
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

    bbp = bp; % base pointer from upper stack
    % predict future values
    Response(:,t+1) = G1*Response(:,t) + C;
    
    EV = Response(bbp:bbp-1+EpsilonPts*ZetaPts*APts,t+1);   
    bbp = bbp+3*EpsilonPts*ZetaPts*APts;
    EV = reshape(EV,EpsilonPts,ZetaPts,APts);
    
    %{
    r(t) = Response(bbp,t);
    bbp = bbp+1;
    w(t) = Response(bbp,t);
    bbp = bbp+1;
    %}
    
    % find KLRatio that clears market
    %{
    options = optimset('Display','iter','TolX',eps);
    KLRatio = fzero(@getKLRatioMetric,[KLRatioMin,KLRatioMax],options);
    %}
    options = optimset('Display','iter','DiffMinChange',1e-3);
    Resi = 1;
    init = KLRatioLast;
    while (Resi>1e-2)
        init = min(max(init+Params.KLRatioBar*randn,KLRatioMin),KLRatioMax);
        [KLRatio,Resi] = fsolve(@getKLRatioMetric,init,options);
    end
    %{
    KLRatio = CoDoSol(Params.KLRatioBar,@getKLRatioMetric,KLRatioMin,KLRatioMax,[1e-8 1e-8]);
    %}
    
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
    AGini(t) = SmltRslt.AGini;
    EntrePopShare(t) = SmltRslt.EntrePopShare;
    KLRatioLast = K(t)/N(t);
    Kp(t) = Kd(t)+K(t);
    ZImplied(t) = Y(t)/(Kp(t)^Gamma*Ns(t)^(1-Gamma));
    DistT{t} = Dist;
    OccPolicyT{t} = VfiRslt.OccPolicy;
    SW(t) = SmltRslt.SW;
    SE(t) = SmltRslt.SE;
    MeanSW(t) = SmltRslt.MeanSW;
    MeanSE(t) = SmltRslt.MeanSE;
    Profit(t) = SmltRslt.Profit;
    MeanProfit(t) = SmltRslt.MeanProfit;
    
    
    %{
    % updates current jump variables
    bbp = bp;
    
    Response(bbp:bbp-1+EpsilonPts*ZetaPts*APts,t) = VfiRslt.EV(:);
    bbp = bbp+EpsilonPts*ZetaPts*APts;
    
    Response(bbp:bbp-1+EpsilonPts*ZetaPts*APts,t) = VfiRslt.VWorker(:);
    bbp = bbp+EpsilonPts*ZetaPts*APts;
    
    Response(bbp:bbp-1+EpsilonPts*ZetaPts*APts,t) = VfiRslt.VManager(:);
    bbp = bbp+EpsilonPts*ZetaPts*APts;
    
    Response(bbp,t) = SmltRslt.r;
    bbp = bbp+1;
    
    Response(bbp,t) = SmltRslt.w;
    bbp = bbp+1;
    %}
    
    %}
    % evolution of jumps
    Response(:,t+1) = Response(:,t+1) + impact*shock(:,t);
    
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
SimulateResult = v2struct(z,G,Lambda,Tax,B,Ks,Kd,K,Ns,Nd,N,r,w,YE,YF,Y,AGini,EntrePopShare,Kp,ZImplied,DistT,OccPolicyT,...
    SW,SE,MeanSW,MeanSE,Profit,MeanProfit);
end