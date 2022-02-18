function vPDC = odePDC(tspan,y0,p,f,pHdata,setup)
% disp(tspan);
% if tspan >= 360
%     disp('uno');
% end

% select initial points
PYR = y0(1); % mM
CO2 = y0(2); % mM
AcAld = y0(3); % mM
ETOH = y0(4); % mM
NADH = y0(5); % mM
NAD = y0(6); % mM

% calulate v (rateEquations)
v_PDC = (p.PDC1_vm.*(PYR./p.PDC1_Kpyr).^p.PDC1_hill)./...
        (1+(PYR./p.PDC1_Kpyr).^p.PDC1_hill);
v_ADH = p.ADH_Vm .* (AcAld .* NADH - ETOH .* NAD ./ p.ADH_Keq);

% mass balances (assumed ideally mixed batch system configuration)
vPDC(1) = - v_PDC;
vPDC(2) = + v_PDC;
vPDC(3) = + v_PDC - v_ADH;
vPDC(4) = + v_ADH;
vPDC(5) = - v_ADH;
vPDC(6) = + v_ADH;

vPDC = vPDC';
end