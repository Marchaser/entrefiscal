% clear;
Params = SETUP;
%{
Params = COMMON(Params);
%}
% CALIBRATE_BETA(Params);
% CALIBRATE(Params);
% LOCAL(Params);
% sys = load('sys');
%}

%{
shock = zeros(3,50);
% shock(3,9) = -0.5;
shock(2,10) = 0.05;
SIMULATE(sys,30,shock,Params);
%}