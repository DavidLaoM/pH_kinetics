function [v] = odeGAPDH_vHeerden(tspan,y0,p,f,data,setup)
%ODEGAPDHR Summary of this function goes here
%   Kinetics from vanHeerden 2014
%   Order of events
%       Select initial concentraitons
%       Calculate v
%       Mass balances
ode_pH = setup.ode_pH;
% setup.ode_pH = 'on';

% select initial points
P3G = y0(1);
ATP = y0(2);
BPG = y0(3);
ADP = y0(4);
NAD = y0(5);
GAP = y0(6);
PHOS = y0(7);
NADH = y0(8);

% calulate v (rateEquations)
switch ode_pH
    case 'on'
        H = 10^(setup.pH_vals(6) - setup.pH_vals(data.i));
        v_GAPDH = (-(p.TDH1_Vmr .* BPG .* NADH .* H ./ (p.TDH1_Kbpg .* p.TDH1_Knadh))...
            + p.TDH1_Vmf .* NAD .* GAP ./ ( p.TDH1_Kgap .* p.TDH1_Knad)) ./...
            ((1 + NAD ./ p.TDH1_Knad + NADH ./ p.TDH1_Knadh) .* (1 + BPG ./...
            p.TDH1_Kbpg + GAP ./ p.TDH1_Kgap));
    otherwise
        v_GAPDH = (-(p.TDH1_Vmr .* BPG .* NADH ./ (p.TDH1_Kbpg .* p.TDH1_Knadh)) + p.TDH1_Vmf .* NAD .* GAP ./ ( p.TDH1_Kgap .* p.TDH1_Knad)) ./ ((1 + NAD ./ p.TDH1_Knad + NADH ./ p.TDH1_Knadh) .* (1 + BPG ./ p.TDH1_Kbpg + GAP ./ p.TDH1_Kgap));
end        
v_PGK = p.PGK_Vm .* (BPG .* ADP - P3G .* ATP ./ p.PGK_Keq);
% v_PGK = p.PGK_Vm .* (BPG .* ADP - P3G .* ATP ./ p.PGK_Keq) ./ (p.PGK_Katp .* p.PGK_Kp3g .* (1 + BPG ./ p.PGK_Kbpg + P3G ./ p.PGK_Kp3g) .* (1 + ADP ./ p.PGK_Kadp + ATP ./ p.PGK_Katp));

% mass balances (assumed ideally mixed batch system configuration)
v(1) = + v_PGK; %P3G
v(2) = + v_PGK; %ATP
v(3) = - v_PGK + v_GAPDH; %BPG
v(4) = - v_PGK; %ADP
v(5) = - v_GAPDH; %NAD
v(6) = - v_GAPDH; %GAP
v(7) = - v_GAPDH; %PHOS
v(8) = + v_GAPDH; % NADH
v=v';

end
%
