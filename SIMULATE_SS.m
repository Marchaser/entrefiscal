function SmltRslt = SIMULATE_SS(z,bond,gov,VfiRslt,Params,Distl,IterMax)
v2struct(Params);
v2struct(VfiRslt);

%{
FullPp = tensor_pchip({AGrid}, cat(1,...
    reshape(ApWorker, 1, EpsilonPts, ZetaPts, APts), ...
    reshape(ApManager, 1, EpsilonPts, ZetaPts, APts)));
FullPp = myppual(FullPp);
%}
FullData = cat(1,...
    reshape(ApWorker, 1, EpsilonPts, ZetaPts, APts), ...
    reshape(ApManager, 1, EpsilonPts, ZetaPts, APts));
FullData = reshape(FullData, [], APts);
FullPp = tensor_pchip({AGrid}, FullData);
FullPp = myppual(FullPp);
%{
pp = struct('form','MKLpp','breaks',{AGrid},...
    'Values',FullData,'coefs',[],'order',[2],...
    'Method',[],'ExtrapolationOrder',[],'thread',8,...
    'orient','curvefit');
FullPp = myppual(pp);
%}

if isequal(Distl,[])
    Distl = 1/EpsilonPts/ZetaPts/APts * ones(EpsilonPts, ZetaPts, APts);
end

ApInterp = myppual(FullPp,AGrid(:)',[],[]);
% ApInterp = myppualMKL_CMEX(int32(NumOfThreads), {AGrid}, FullPp.coefs, [], int32([2]), int32(2*EpsilonPts*ZetaPts), [], [AGrid(:)'], [], [], []);
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

Err = 1;
count = 0;
if isequal(IterMax,[])
    IterMax = inf;
end
Dist = Distl;
while (Err>TolEqSs && count<IterMax)
    count = count+1;
    
    Distl = Dist;
    Dist = update_dist(Distl,OccPolicy,ApManagerCellInt,ApWorkerCellInt,ApManagerLeftShare,ApWorkerLeftShare,EpsilonTrans,ZetaTrans);

    Err = sum(abs(Dist(:)-Distl(:)));
    
    if Params.ShowDetail==1
        if mod(count,100)==1
            display(['Iter: ' num2str(count) ', Err: ' num2str(Err)]);
        end
    end
end

Ks = sum(Distl(:).*AMesh(:));
Kd = sum(Distl(:).*OccPolicy(:).*KManager(:));
K = Ks - Kd - bond;
Nd = sum(Distl(:).*OccPolicy(:).*NManager(:));
Ns = sum(Distl(:).*(1-OccPolicy(:)).*NWorker(:).*EpsilonMesh(:));
N = Ns - Nd;
KLRatio = K/N;
r = z*Gamma*KLRatio^(Gamma-1)-Delta;
w = z*(1-Gamma)*KLRatio^Gamma;
ADist = sum(reshape(Distl,[],APts));
AGini = gini_dist(ADist,AGrid);
EntrePopShare = sum(Distl(:).*OccPolicy(:));
WorkerPopShare = 1-EntrePopShare;
MeanHours = sum(Distl(:).*(1-OccPolicy(:)).*NWorker(:))/WorkerPopShare;
MeanNs = Ns/WorkerPopShare;
MeanFirmSize = Nd/EntrePopShare/MeanNs;
KShare = Kd/Ks;
NShare = Nd/Ns;

% employment distribution
ManagerHasEmp = sum(Distl(:).*(NManager(:)/MeanNs>1).*OccPolicy(:));
EmpSize1To5 = sum(Distl(:).*(NManager(:)/MeanNs<=5 & NManager(:)/MeanNs>1).*OccPolicy(:))/ManagerHasEmp;
EmpSize5To10 = sum(Distl(:).*(NManager(:)/MeanNs<=10 & NManager(:)/MeanNs>5).*OccPolicy(:))/ManagerHasEmp;
EmpSize10To20 = sum(Distl(:).*(NManager(:)/MeanNs<=20 & NManager(:)/MeanNs>10).*OccPolicy(:))/ManagerHasEmp;
EmpSize20To100 = sum(Distl(:).*(NManager(:)/MeanNs<=100 & NManager(:)/MeanNs>20).*OccPolicy(:))/ManagerHasEmp;
EmpSizeAbove100 = sum(Distl(:).*(NManager(:)/MeanNs>100).*OccPolicy(:))/ManagerHasEmp;

CE = sum(Distl(:).*CManager(:));
CW = sum(Distl(:).*CWorker(:));
C = CE+CW;
YE = sum(Distl(:).*YManager(:).*OccPolicy(:));
YF = z*K^Gamma*N^(1-Gamma);
Y = YE+YF;
KYRatio = Ks/Y;

% investment
IE = Delta*Kd;
IF = Delta*K;
I = IE+IF;

Tax = Rho0+Rho1*bond+Rho2*gov;
Bp = bond*(1+r)+gov-Tax;

% exit rate
%{
Survival = zeros(size(Dist));
for iE=1:EpsilonPts
    for iZeta=1:ZetaPts
        for iA=1:APts
            iApManagerCell = ApManagerCell(iE,iZeta,iA);
            iApManagerLeftShare = ApManagerLeftShare(iE,iZeta,iA);
            for iEPrime=1:EpsilonPts
                for iZetaPrime=1:ZetaPts
                    Survival(iE,iZeta,iA) = Survival(iE,iZeta,iA) + ...
                        (OccPolicy(iEPrime,iZetaPrime,iApManagerCell)*iApManagerLeftShare ...
                        +OccPolicy(iEPrime,iZetaPrime,iApManagerCell+1)*(1-iApManagerLeftShare)) ...
                        *EpsilonTrans(iE,iEPrime)*ZetaTrans(iZeta,iZetaPrime);
                end
            end
        end
    end
end
%}
SurvivalPolicy = trans_occ(OccPolicy,ApManagerCellInt,ApManagerLeftShare,EpsilonTrans,ZetaTrans);
Survival = sum(Distl(:).*SurvivalPolicy(:).*OccPolicy(:))/EntrePopShare;

SmltRslt = v2struct(Dist,Ks,K,Kd,Nd,Ns,N,KLRatio,r,w,Tax,Bp,YE,YF,Y,IE,IF,I,EntrePopShare,WorkerPopShare,MeanHours,MeanFirmSize,...
    AGini,KShare,NShare,Survival,KYRatio,EmpSize1To5,EmpSize5To10,EmpSize10To20,EmpSize20To100,EmpSizeAbove100);
end
