function [vObs, vextra] = calcRates_gapdhFWDREV(xtemp,simResult,data,setup)
enzymeName = setup.enzymeName;
constantVm = setup.constantVm;
ode = setup.ode;
sourceVm = setup.sourceVm;

% metabolites
P3G = simResult.y(:,1);
ATP = simResult.y(:,2);
BPG = simResult.y(:,3);
ADP = simResult.y(:,4);
NAD = simResult.y(:,5);
GAP = simResult.y(:,6);
PHOS = simResult.y(:,7);
NADH = simResult.y(:,8);

% recall parameters
p.TDH1_Kgap=10.^xtemp(1).*2.48; % mM
p.TDH1_Kbpg=10.^xtemp(2).*1.18; % mM
p.TDH1_Knad=10.^xtemp(3).*2.92; %mM
p.TDH1_Knadh=10.^xtemp(4).*0.022; % mM
p.TDH1_Vmf = data.chosenVmaxFWD * 10 .^ xtemp(5); %(UNIT!) 
p.TDH1_Vmr = data.chosenVmaxREV * 10 .^ xtemp(6); %(UNIT!) 
p.GAPDH_Keq = data.chosenKeqGAPDH;
p.PGK_Vm = 1306.45 / 60;
p.PGK_Keq = data.chosenKeqPGK;      

% calculate the rates (adjusted)
v_GAPDH = (p.TDH1_Vmf .* GAP .* NAD ./ (p.TDH1_Kgap .* p.TDH1_Knad) - ...
    p.TDH1_Vmr .* BPG .* NADH ./ (p.TDH1_Kbpg .* p.TDH1_Knadh)) ./ ...
    ((1 + NAD ./ p.TDH1_Knad + NADH ./ p.TDH1_Knadh) .* (1 + BPG ./...
    p.TDH1_Kbpg + GAP ./ p.TDH1_Kgap));
vObs = abs(v_GAPDH);
vextra = 1;

end