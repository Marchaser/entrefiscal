function [VfiRslt, Flag] = VFI_SS(z,r,w,tax,lambda,tauL,tauR,tauPi,Params,EV,IterMax)
v2struct(Params);

% worker
Worker.data = repmat([Sigma;Chi;tauL;tauR;r;w;tax],1,NX);
WorkerBudget = AMesh*(1+r*(1-tauR)) + w*EpsilonMesh*(1-tauL) - tax;
Worker.data = [
    Worker.data;
    EpsilonMesh(:)';
    AMesh(:)';
    WorkerBudget(:)';
    ];
Worker.ub = min(max(WorkerBudget(:)',Worker.lb),AMax);

% manager
A = z*ZetaMesh;
Coef1 = (1-Gamma)*(r+Delta) / (Gamma*w);
KManager = power(A*Theta*Gamma*power(Coef1,Theta*(1-Gamma))/(r+Delta), 1/(1-Theta));
KManager = min(KManager,AMesh*lambda);
NManager = Coef1*KManager;
YManager = A.*KManager.^(Gamma*Theta).*NManager.^((1-Gamma)*Theta);
PiManager = YManager - r*(KManager-AMesh) - Delta*KManager - w*NManager;
ManagerBudget = PiManager*(1-tauPi) + AMesh - tax;
Manager.data = repmat([Sigma Chi Theta Gamma lambda Delta tauPi tauR z r w tax]',[1 NX]);
Manager.data = [
    Manager.data;
    ZetaMesh(:)';
    AMesh(:)';
    ManagerBudget(:)';
    ];
Manager.ub = min(max(ManagerBudget(:)',Manager.lb),AMax);

% value
if isequal(EV,[])
    EV = zeros(EpsilonPts, ZetaPts, APts);
end
% Worker.x0 = (Worker.lb + Worker.ub)/2;
% Manager.x0 = (Manager.lb + Manager.ub)/2;
Worker.x0 = Worker.lb+1e-2;
Manager.x0 = Manager.lb+1e-2;

%
if isequal(IterMax,[])
    IterMax = inf;
end
Err = 1;
count = 0;
while(Err > TolVfi && count < IterMax)
    count = count+1;
    % interpolation
    APp = tensor_pchip(AGrid, reshape(EV, [], APts));
    AMklPp = myppual(APp);
    %}
    %{
    AMklPp.breaks = {AGrid};
    AMklPp.order = int32(2);
    AMklPp.dim = int32(EpsilonPts*ZetaPts);
    AMklPp.coefs = myppualMKL_CMEX(int32(-NumOfThreads), {AGrid}, permute(reshape(EV, [], APts),[2 1]), [], int32([2]), int32(EpsilonPts*ZetaPts), [], [], [], [], []);
    %}
    
    % worker's problem
    Worker.pp = AMklPp;
    [XWorker,VWorker,EFWorker] = dpopt_worker(Worker);
    PolicyWorker = reshape(Worker.rhs_return(1:3,:), [3 EpsilonPts ZetaPts APts]);
    
    % manager's problem
    Manager.pp = AMklPp;
    [XManager,VManager,EFManager] = dpopt_manager(Manager);
    PolicyManager = reshape(Manager.rhs_return(1:3,:), [3 EpsilonPts ZetaPts APts]);
    CManager = PolicyManager(2,:);
    
    VManager = reshape(VManager,[EpsilonPts ZetaPts APts]);
    VWorker = reshape(VWorker,[EpsilonPts ZetaPts APts]);
    
    %{
    V = VManager + normcdf((VWorker-VManager)/USigma).*(VWorker-VManager) ...
        +USigma*normpdf((VWorker-VManager)/USigma);
    V(:,1,:) = VWorker(:,1,:);
    %}
    V = VWorker + log(1 + exp(VManager-VWorker));
    V(:,1,:) = VWorker(:,1,:);
    
    EVNew = Beta*EpsilonZetaTrans * reshape(V, EpsilonPts*ZetaPts, APts);
    Err = max(abs(EV(:)-EVNew(:)));
    EV = reshape(EVNew, [EpsilonPts ZetaPts APts]);
    
    if Params.ShowDetail==1
        if mod(count,100)==1
            display(['Iter: ' num2str(count) ', Err: ' num2str(Err)]);
        end
    end
    %}
end
% OccPolicy = exp(VManager) ./ (exp(VManager)+exp(VWorker));
% VManagerSmall = min(VManager(CManager>0));
% VManager(CManager<0) = VManagerSmall-100;
% OccPolicy = normcdf((VManager-VWorker)/USigma);
OccPolicy = 1 ./ (1+exp(VWorker-VManager));
OccPolicy = reshape(OccPolicy, EpsilonPts, ZetaPts, APts);
% OccPolicy(CManager<0) = 0;
% CManager(CManager<0) = 0;
OccPolicy(:,1,:) = 0;

ApWorker = PolicyWorker(1,:);
NWorker = PolicyWorker(2,:);
CWorker = PolicyWorker(3,:);
ApManager = PolicyManager(1,:);

VfiRslt = v2struct(EV,V,VWorker,VManager,KManager,NManager,YManager,PiManager,ApWorker,NWorker,CWorker,...
    ApManager,CManager,OccPolicy);
Flag = 1;
end

