function vPGI = odePGI(tspan,y0,p,f,pHdata,setup)
% disp(tspan);
% if tspan >= 360
%     disp('uno');
% end

% select initial points
G6P = y0(1); % mM
F6P = y0(2); % mM
ATP = y0(3); % mM
FBP = y0(4); % mM
ADP = y0(5); % mM
DHAP = y0(6); % mM
GAP = y0(7); % mM
NADH = y0(8); % mM
G3P = y0(9); % mM
NAD = y0(10); % mM

% calulate v (rateEquations)
v_PGI = ((p.PGI1_Vm).*(G6P-(F6P./p.PGI1_Keq)))./...
    (1+G6P./p.PGI1_Kg6p+1+F6P./p.PGI1_Kf6p-1);
v_PFK = p.PFK1_Vm .* (F6P .* ATP - FBP .* ADP ./ p.PGI1_Keq);
v_ALD = p.FBA1_Vm .* (FBP - (GAP .* DHAP) ./ p.FBA1_Keq);
v_TPI = p.TPI1_Vm .* (DHAP - GAP ./ p.TPI1_Keq);
v_GPD = p.GPD1_Vm .* (DHAP .* NADH - G3P .* NAD ./ p.GPD1_Keq);

% mass balances (assumed ideally mixed batch system configuration)
vPGI(1) = - v_PGI;
vPGI(2) = + v_PGI - v_PFK;
vPGI(3) = - v_PFK;
vPGI(4) = + v_PFK - v_ALD;
vPGI(5) = + v_PFK;
vPGI(6) = + v_ALD - v_TPI - v_GPD;
vPGI(7) = + v_ALD + v_TPI;
vPGI(8) = - v_GPD;
vPGI(9) = + v_GPD;
vPGI(10) = + v_GPD;

vPGI = vPGI';
end