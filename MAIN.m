% clear;
Params = SETUP;
%{
Params = COMMON(Params);
%}
% CALIBRATE_BETA(Params);
% LOCAL(Params);
% sys = load('sys');

shock = zeros(3,50);
shock(2,1) = 0.05;
SIMULATE(sys,20,shock,Params);