function vTPI = odeTPI(tspan,y0,p,f,pHdata,setup)
% disp(tspan);
% if tspan >= 360
%     disp('uno');
% end

% select initial points
FBP = y0(1); % mM % just write down zero there
DHAP = y0(2); % mM
GAP = y0(3); % mM
G3P = y0(4); % mM
NADH = y0(5); % mM
NAD = y0(6); % mM

% calulate v (rateEquations)
% v_ALD = p.FBA1_Vm .* (FBP-(GAP.*DHAP)./p.FBA1_Keq);
v_ALD = 0;
v_GPD = p.GPD1_Vm .* (DHAP .* NADH - G3P .* NAD ./ p.GPD1_Keq);
v_TPI = (p.TPI1_Vm./p.TPI1_Kdhap.*(DHAP-GAP./p.TPI1_Keq))./...
    (1+DHAP./p.TPI1_Kdhap+GAP./p.TPI1_Kglyceral3p);
% v_TPI = p.TPI1_Vm .* (0 - 1) * 2;
% v_TPI = p.TPI1_Vm .* (0 - 1);

% mass balances (assumed ideally mixed batch system configuration)
vTPI(1) = - v_ALD;
vTPI(2) = + v_ALD - v_GPD - v_TPI;
vTPI(3) = + v_ALD + v_TPI;
vTPI(4) = + v_GPD;
vTPI(5) = - v_GPD;
vTPI(6) = + v_GPD;
% % % % vTPI(1) = - v_ALD;
% % % % vTPI(2) = + v_ALD - v_GPD - 2* v_TPI;
% % % % vTPI(3) = + v_ALD + 2* v_TPI;
% % % % vTPI(4) = + v_GPD;
% % % % vTPI(5) = - v_GPD;
% % % % vTPI(6) = + v_GPD;

vTPI=vTPI';
end