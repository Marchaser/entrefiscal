function CALIBRATE(Params)
x0 = Params.CaliX0;
Params.TolOpt = 1e-10;
Params.TolVfi = 1e-8;
Params.TolEqSs = 1e-8;
Params.ShowDetail = 1;
%{
x0(1) = Params.Beta;
x0(2) = Params.Delta;
x0(3) = Params.Chi;
x0(4) = Params.EntreShock1;
x0(5) = Params.EntreShock2;
x0(6) = Params.USigma;
x0(7) = Params.ZetaMu;
x0(8) = Params.ZetaSigma;
x0(9) = Params.Theta;
x0(10) = Params.Lambda;
x0(11) = Params.BBar;
x0(12) = Params.GBar;
%}

LastX = x0';
bound = [
    0 1
    0 0.1
    0 1
    0 1
    0 1
    -inf inf
    1e-6 1
    0 1
    1 inf
    0 inf
    0 inf
    ];
XLb = bound(:,1)';
XUb = bound(:,2)';

% load('VfiRslt');
% EV = VfiRslt.EV;
EV = [];
Dist = [];

    function MomentsMetric = ComputeMomentsMetric(x)
        display('Current Parameters:');
        display(x');
        IsComputingDiff = sum(x~=LastX)==1;
        
        NewParams = Params;
        
        NewParams.Beta = x(1);
        NewParams.Delta = x(2);
        NewParams.Chi = x(3);
        NewParams.EntreShock1 = x(4);
        NewParams.EntreShock2 = x(5);
        NewParams.ZetaMu = x(6);
        NewParams.ZetaSigma = x(7);
        NewParams.Theta = x(8);
        NewParams.LambdaBar = x(9);
        BBar = x(10);
        GBar = x(11);
        Tax = BBar*Params.RBar + GBar;
        
        NewParams = COMMON(NewParams);
        DataMoments = [
            Params.RBar
            10.6
            0.33
            0.0755
            0.81
            0.7568
            0.146
            0.178
            0.218
            0.2086
            0.3000*4
            ]';
        
        R = Params.RBar;
        KLRatio = ((R+NewParams.Delta)/Params.Gamma)^(1/(Params.Gamma-1));
        W = (1-Params.Gamma)*KLRatio^Params.Gamma;
        display('VFI...');
        VfiRslt = VFI_SS(NewParams.ZBar,R,W, ...
            Tax,NewParams.LambdaBar,NewParams.TauLBar,NewParams.TauRBar,NewParams.TauPiBar,NewParams,EV,[]);
        if ~IsComputingDiff
            EV = VfiRslt.EV;
            save('VfiRslt', 'VfiRslt');
            LastX = x;
        end
        
        display('Simulate...');
        SmltRslt = SIMULATE_SS(NewParams.ZBar,BBar,GBar,VfiRslt,NewParams,Dist,[]);
        if ~IsComputingDiff
            Dist = SmltRslt.Dist;
        end
        
        ModelMoments = [
            SmltRslt.r
            SmltRslt.KYRatio
            SmltRslt.MeanHours
            SmltRslt.EntrePopShare
            SmltRslt.AGini
            SmltRslt.EmpSize1To5
            SmltRslt.EmpSize5To10
            SmltRslt.KShare
            SmltRslt.NShare
            Tax / SmltRslt.Y
            BBar / SmltRslt.Y
            ]';
            
        display('Current Moments: ');
        display([ModelMoments;DataMoments]);
        MomentsMetric = (ModelMoments - DataMoments)./DataMoments;
        %{
        MomentsMetric(1) = MomentsMetric(1)*10;
        MomentsMetric(4) = MomentsMetric(4)*10;
        MomentsMetric(end-3:end-2) = MomentsMetric(end-3:end-2)*10;
        %}
    end

% options = [];
% LastX = CoDoSol(x0', @(x) ComputeMomentsMetric(x)', XLb', XUb', [1e-6, 1e-6], options)';
% fsolve(@(x) ComputeMomentsMetric(x)', x0', optimoptions('fsolve','Display','Iter','FinDiffRelStep',1e-3,'DiffMinChange',1e-6,'TypicalX',x0'));

options = optimoptions('lsqnonlin','Display','Iter','FinDiffRelStep',1e-3,'DiffMinChange',1e-6,'TypicalX',x0');
lsqnonlin(@(x)ComputeMomentsMetric(x)', x0', XLb', XUb', options);

end