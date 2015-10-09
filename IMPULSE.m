Params = SETUP;
sys = load('sys');

% government sending shock, starting from period 2
shock = zeros(3,50);
shock(2,2) = Params.GBar*0.1;
GImpulse = SIMULATE(sys,30,shock,Params);

% credit tight shock, starting from period 1
shock = zeros(3,50);
shock(3,1) = 0.26*Params.LambdaBar;
LambdaImpulse = SIMULATE(sys,30,shock,Params);

% credit tight shock, following government spending shock
shock = zeros(3,50);
shock(3,1) = 0.26*Params.LambdaBar;
shock(2,2) = Params.GBar*0.1;
GLambdaImpulse = SIMULATE(sys,30,shock,Params);

shock = zeros(3,50);
shock(3,1) = -0.26*Params.LambdaBar;
MinusLambdaImpulse = SIMULATE(sys,30,shock,Params);

shock = zeros(3,50);
shock(3,1) = -0.26*Params.LambdaBar;
shock(2,2) = Params.GBar*0.1;
GMinusLambdaImpulse = SIMULATE(sys,30,shock,Params);

% plot distribution
Params = SETUP;
Params = COMMON(Params);
[ZetaMesh,AMesh] = ndgrid(Params.ZetaGrid,Params.AGrid);

%


% state dependence?
figure;
hold on;
plot(GImpulse.EntrePopShare-GImpulse.EntrePopShare(1),'k');
plot(GLambdaImpulse.EntrePopShare-LambdaImpulse.EntrePopShare);
hold off;

figure;
hold on;
plot(GImpulse.Y-GImpulse.Y(1),'k');
plot(GLambdaImpulse.Y-LambdaImpulse.Y);
hold off;
