function vALD = odeALD(tspan,y0,p,f,pHdata,setup)
% disp(tspan);
% if tspan >= 360
%     disp('uno');
% end

% select initial points
FBP = y0(1); % mM
DHAP = y0(2); % mM
GAP = y0(3); % mM
G3P = y0(4); % mM
NADH = y0(5); % mM
NAD = y0(6); % mM

% calulate v (rateEquations)
% v_ALD = p.FBA1_Vm .* (FBP-(GAP.*DHAP)./p.FBA1_Keq);
v_ALD = p.FBA1_Vm .* (FBP-(GAP.*DHAP)./p.FBA1_Keq) ./...
    (p.FBA1_Kf16bp .* (FBP./p.FBA1_Kf16bp + (1 + GAP./p.FBA1_Kglyceral3p).*(1 + DHAP./p.FBA1_Kdhap)));
v_GPD = p.GPD1_Vm .* (DHAP .* NADH - G3P .* NAD ./ p.GPD1_Keq);
v_TPI = p.TPI1_Vm .* (DHAP - GAP ./ p.TPI1_Keq);
%test
% v_ALD = p.FBA1_Vm;
% v_GPD = 100000 .* (DHAP .* NADH - G3P .* NAD ./ p.GPD1_Keq);
% v_TPI = 100000 .* (DHAP - GAP ./ p.TPI1_Keq);

% mass balances (assumed ideally mixed batch system configuration)
vALD(1) = - v_ALD;
vALD(2) = + v_ALD - v_GPD - v_TPI;
vALD(3) = + v_ALD + v_TPI;
vALD(4) = + v_GPD;
vALD(5) = - v_GPD;
vALD(6) = + v_GPD;

vALD=vALD';
end