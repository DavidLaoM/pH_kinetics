function vPGM = odePGM(tspan,y0,p,f,pHdata,setup)
% disp(tspan);
% if tspan >= 360
%     disp('uno');
% end

% select initial points
P3G = y0(1); % mM
P2G = y0(2); % mM
PEP = y0(3); % mM
PYR = y0(4); % mM
LAC = y0(5); % mM
ADP = y0(6); % mM
ATP = y0(7); % mM
NADH = y0(8); % mM
NAD = y0(9); % mM

% calulate v (rateEquations)
v_PGM = ((p.GPM1_vm ./ p.GPM1_K3pg) .* (P3G - P2G ./ p.GPM1_Keq))./...
    (1 + P3G ./ p.GPM1_K3pg + P2G ./ p.GPM1_K2pg);
v_ENO = p.ENO_Vm .* (P2G - PEP ./ p.ENO_Keq);
v_PYK = p.PYK_Vm .* (PEP .* ADP - PYR .* ATP ./ p.PYK_Keq);
v_LDH = p.LDH1_Vm .*(PYR .* NADH - LAC .* NAD ./ p.LDH1_Keq);

% mass balances (assumed ideally mixed batch system configuration)
vPGM(1) = - v_PGM;
vPGM(2) = + v_PGM - v_ENO;
vPGM(3) = + v_ENO - v_PYK;
vPGM(4) = + v_PYK - v_LDH;
vPGM(5) = + v_LDH;
vPGM(6) = - v_PYK;
vPGM(7) = + v_PYK;
vPGM(8) = - v_LDH;
vPGM(9) = + v_LDH;

vPGM = vPGM';
end