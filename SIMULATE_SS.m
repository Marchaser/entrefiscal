function SmltRslt = SIMULATE_SS(bond,gov,VfiRslt,Params,Dist,IterMax)
v2struct(Params);
v2struct(VfiRslt);

FullPp = tensor_pchip({AGrid}, cat(1,...
    reshape(ApWorker, 1, EpsilonPts, ZetaPts, APts), ...
    reshape(ApManager, 1, EpsilonPts, ZetaPts, APts)));
FullPp = myppual(FullPp);

if isequal(Dist,[])
    Dist = 1/EpsilonPts/ZetaPts/APts * ones(EpsilonPts, ZetaPts, APts);
end

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

Err = 1;
count = 0;
while (Err>TolEqSs)
    count = count+1;

    
    TDist = update_dist(Dist,OccPolicy,ApManagerCellInt,ApWorkerCellInt,ApManagerLeftShare,ApWorkerLeftShare,EpsilonTrans,ZetaTrans);
    
    % update distribution
    % prototype
    Err = sum(abs(TDist(:)-Dist(:)));
    Dist = TDist;
    
    if mod(count,100)==1
        display(['Iter: ' num2str(count) ', Err: ' num2str(Err)]);
    end
    
    if (Err<=TolEqSs)
    end
end

ADist = sum(reshape(Dist,[],APts));
AGini = gini_dist(ADist,AGrid);
Ks = sum(Dist(:).*AMesh(:));
Kd = sum(Dist(:).*OccPolicy(:).*KManager(:));
K = Ks - Kd - bond;
Nd = sum(Dist(:).*OccPolicy(:).*NManager(:));
Ns = sum(Dist(:).*(1-OccPolicy(:)).*NWorker(:).*EpsilonMesh(:));
N = Ns - Nd;
KLRatio = K/N;
r = Gamma*KLRatio^(Gamma-1)-Delta;
w = (1-Gamma)*KLRatio^Gamma;
EntrePopShare = sum(Dist(:).*OccPolicy(:));
WorkerPopShare = 1-EntrePopShare;
MeanHours = sum(Dist(:).*(1-OccPolicy(:)).*NWorker(:))/WorkerPopShare;
MeanNs = Ns/WorkerPopShare;
MeanFirmSize = Nd/EntrePopShare/MeanNs;
KShare = Kd/Ks;
NShare = Nd/Ns;

% employment distribution
ManagerHasEmp = sum(Dist(:).*(NManager(:)/MeanNs>1).*OccPolicy(:));
EmpSize1To5 = sum(Dist(:).*(NManager(:)/MeanNs<=5 & NManager(:)/MeanNs>1).*OccPolicy(:))/ManagerHasEmp;
EmpSize5To10 = sum(Dist(:).*(NManager(:)/MeanNs<=10 & NManager(:)/MeanNs>5).*OccPolicy(:))/ManagerHasEmp;
EmpSize10To20 = sum(Dist(:).*(NManager(:)/MeanNs<=20 & NManager(:)/MeanNs>10).*OccPolicy(:))/ManagerHasEmp;
EmpSize20To100 = sum(Dist(:).*(NManager(:)/MeanNs<=100 & NManager(:)/MeanNs>20).*OccPolicy(:))/ManagerHasEmp;
EmpSizeAbove100 = sum(Dist(:).*(NManager(:)/MeanNs>100).*OccPolicy(:))/ManagerHasEmp;

CE = sum(Dist(:).*CManager(:));
CW = sum(Dist(:).*CWorker(:));
C = CE+CW;
YE = sum(Dist(:).*YManager(:).*OccPolicy(:));
YF = K^Gamma*N^(1-Gamma);
Y = YE+YF;
KYRatio = Ks/Y;

% investment
IE = Delta*Kd;
IF = Delta*K;
I = IE+IF;

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
Survival = sum(Dist(:).*SurvivalPolicy(:).*OccPolicy(:))/EntrePopShare;
SmltRslt = v2struct(Dist,Ks,K,Kd,Nd,Ns,N,KLRatio,r,w,YE,YF,Y,IE,IF,I,EntrePopShare,WorkerPopShare,MeanHours,MeanFirmSize,...
    AGini,KShare,NShare,Survival,KYRatio,EmpSize1To5,EmpSize5To10,EmpSize10To20,EmpSize20To100,EmpSizeAbove100);
end
