function CALIBRATE_BETA(Params)
Params = COMMON(Params);
x0 = Params.Beta;
LastX = x0';
bound = [
    0 1
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
        DataMoments = [
            Params.RBar
            ]';
        Tax = Params.BBar*Params.RBar + Params.GBar;
        
        R = Params.RBar;
        KLRatio = ((R+NewParams.Delta)/Params.Gamma)^(1/(Params.Gamma-1));
        W = (1-Params.Gamma)*KLRatio^Params.Gamma;
        display('VFI...');
        VfiRslt = VFI_SS(NewParams.ZBar,R,W, ...
            Tax,NewParams.LambdaBar,NewParams.TauLBar,NewParams.TauRBar,NewParams.TauPiBar,NewParams,EV,[]);
        %         if ~IsComputingDiff
        EV = VfiRslt.EV;
        save('VfiRslt', 'VfiRslt');
        LastX = x;
        %         end
        
        display('Simulate...');
        SmltRslt = SIMULATE_SS(Params.BBar,Params.GBar,VfiRslt,NewParams,Dist,[]);
        if ~IsComputingDiff
            Dist = SmltRslt.Dist;
        end
        
        ModelMoments = [
            SmltRslt.r
            ]';
        
        display('Current Moments: ');
        display([ModelMoments;DataMoments]);
        MomentsMetric = (ModelMoments - DataMoments)./DataMoments;
    end

% options = [];
% LastX = CoDoSol(x0', @(x) ComputeMomentsMetric(x)', XLb', XUb', [1e-8, 1e-8], options)';
% fsolve(@(x) ComputeMomentsMetric(x)', x0', optimoptions('fsolve','Display','Iter','FinDiffRelStep',1e-4,'DiffMinChange',1e-6));

options = optimoptions('lsqnonlin','Display','Iter','FinDiffRelStep',1e-6,'DiffMinChange',1e-6,'TypicalX',x0','TolX',1e-8,'TolFun',1e-20);
lsqnonlin(@(x)ComputeMomentsMetric(x)', x0', XLb', XUb', options);

end