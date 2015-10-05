Params = SETUP;
Params = COMMON(Params);
v2struct(Params);
load('VfiRslt');
% VfiRslt =
% VFI_SS(ZBar,RBar,WBar,TaxBar,TauLBar,TauRBar,TauPiBar,Params,VfiRslt.V,[]);
VfiRslt = VFI_SS(ZBar,RBar,WBar,TaxBar,LambdaBar,TauLBar,TauRBar,TauPiBar,Params,VfiRslt.EV,[]);
save('VfiRslt','VfiRslt');
% load('SmltRslt');
% SmltRslt = SIMULATE_SS(BondBar,GBar,VfiRslt,Params,SmltRslt.Dist,[]);
SmltRslt = SIMULATE_SS(BBar,GBar,VfiRslt,Params,[],[]);
% save('SmltRslt','SmltRslt');

LambdaBar = LambdaBar+1e-5;
VfiRslt2 = VFI_SS(ZBar,RBar,WBar,TaxBar,LambdaBar,TauLBar,TauRBar,TauPiBar,Params,VfiRslt.EV,[]);
SmltRslt2 = SIMULATE_SS(BBar,GBar,VfiRslt2,Params,[],[]);

Params = SETUP;
Params = COMMON(Params);
v2struct(Params);
% RIdx = (EpsilonPts*ZetaPts*APts)*(1+2