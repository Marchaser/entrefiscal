function SIMULATE2(sys,T,shock,Params)
Params = COMMON(Params);
% Params.ShowDetail = 1;

G1 = sys.G1;
C = sys.C;
Ss = sys.Ss(:);
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

Response(:,1) = Ss;

Dist = zeros(EpsilonPts,ZetaPts,APts);
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
    
    Err = 1;
    Tol = 1e-6;
    speed = 0.05;
    T_Response = Response(:,t); % states are always left untouched
    while (Err>Tol)
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
        
        % solve decision problem and aggregates
        % solve decision problem
        VfiRslt = VFI_SS(z(t),r(t),w(t),Tax(t),Lambda(t),TauLBar,TauRBar,TauPiBar,Params,EV,1);
        % simulate
        SmltRslt = SIMULATE_SS(z(t),B(t),G(t),VfiRslt,Params,Dist,1);
        
        % get implied aggregates
        bbp = bp;
        T_Response(bbp:bbp-1+EpsilonPts*ZetaPts*APts) = VfiRslt.VWorker(:);
        bbp = bbp+EpsilonPts*ZetaPts*APts;
        
        T_Response(bbp:bbp-1+EpsilonPts*ZetaPts*APts) = VfiRslt.VManager(:);
        bbp = bbp+EpsilonPts*ZetaPts*APts;
        
        T_Response(bbp) = SmltRslt.r;
        bbp = bbp+1;
        
        T_Response(bbp) = SmltRslt.w;
        bbp = bbp+1;        
  
        Err = abs(max(T_Response(:)-Response(:,t)));
        display(Err);
        
        Response(:,t) = T_Response(:)*speed + Response(:,t)*(1-speed);
    end
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
    EntrePopShare(t) = SmltRslt.EntrePopShare;
    
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