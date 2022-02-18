function [v] = odeGAPDHr(tspan,y0,p,f,data,setup)
%ODEGAPDHR Summary of this function goes here
%   Specific kinetics and best initial guesses is still a bit unclear.

% select initial points
P3G = y0(1);
ATP = y0(2);
BPG = y0(3);
ADP = y0(4);
NAD = y0(5);
GAP = y0(6);
Pi = y0(7);
NADH = y0(8);

% % % extra check for non-negative values: if negative values are reached, stay
% % % at zero.
% % if P3G <= 0, P3G = 0; end
% % if ATP <= 0, ATP = 0; end
% % if BPG <= 0, BPG = 0; end
% % if ADP <= 0, ADP = 0; end
% % if NAD <= 0, NAD = 0; end
% % if GAP <= 0, GAP = 0; end
% % if Pi <= 0, Pi = 0; end
% % if NADH <= 0, NADH = 0; end
% % % % 

% calulate v (rateEquations)
v_PGK = 1000 .* P3G * p.linkreaction;
v_GAPDHr = p.TDH1_Vm .* ( (BPG .* NADH) - (GAP .* NAD) ./ p.TDH1_Keq ) ./ (p.TDH1_Kbpg .* p.TDH1_Knadh) ./ ( (1 + GAP./p.TDH1_Kgap + BPG./p.TDH1_Kbpg).*(1 + NAD./p.TDH1_Knad +NADH./p.TDH1_Knadh) );
% % % % v_GAPDH = p.TDH1_Vm .* (NAD .* GAP - BPG .* NADH ./ p.TDH1_Keq) ./ ((p.TDH1_Kgap .* p.TDH1_Knad) .* (1 + NAD ./ p.TDH1_Knad + NADH ./ p.TDH1_Knadh) .* (1 + BPG ./ p.TDH1_Kbpg + GAP ./ p.TDH1_Kgap));
% % % % v_PGK = p.linkreaction .* (BPG .* ADP - P3G .* ATP ./ 3200);


% mass balances (assumed ideally mixed batch system configuration)
v(1) = - v_PGK;
v(2) = - v_PGK;
v(3) = + v_PGK - v_GAPDHr;
v(4) = + v_PGK;
v(5) = + v_GAPDHr;
v(6) = + v_GAPDHr;
v(7) = + v_GAPDHr;
v(8) = - v_GAPDHr;

% % % % v(1) = + v_PGK; %P3G
% % % % v(2) = + v_PGK; %ATP
% % % % v(3) = - v_PGK + v_GAPDH; %BPG
% % % % v(4) = - v_PGK; %ADP
% % % % v(5) = - v_GAPDH; %NAD
% % % % v(6) = - v_GAPDH; %GAP
% % % % v(7) = - v_GAPDH; %PHOS
% % % % v(8) = + v_GAPDH; % NADH


% % % % extra check for non-negative values: if negative values are reached, stay
% % % % at zero.
% % % if P3G <= 0, v(1) = 0; end
% % if ATP <= 0, v(2) = 0; end
% % % if BPG <= 0, v(3) = 0; end
% % % if ADP <= 0, v(4) = 0; end
% % % if NAD <= 0, v(5) = 0; end
% % % if GAP <= 0, v(6) = 0; end
% % % if Pi <= 0, v(7) = 0; end
% % % if NADH <= 0, v(8) = 0; end
% % % % % 

v=v';
end
%
