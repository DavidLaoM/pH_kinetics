function vPYK = odePYK(tspan,y0,p,f,pHdata,setup)
% disp(tspan);
% if tspan >= 360
%     disp('uno');
% end

% select initial points
ADP = y0(1); % mM
NADH = y0(2); % mM
FBP = y0(3); % mM
PEP = y0(4); % mM
ATP = y0(5); % mM
NAD = y0(6); % mM
PYR = y0(7); % mM
LAC = y0(8); % mM

% calulate v (rateEquations)
v_PYK=(((p.PYK1_Vm./...
    (p.PYK1_Kadp.*p.PYK1_Kpep).*ADP.*PEP)./((1+ADP./p.PYK1_Kadp).*...
    (1+PEP./p.PYK1_Kpep))).*((PEP./p.PYK1_Kpep+1).^p.PYK1_hill./...
    (p.PYK1_L.*((ATP./p.PYK1_Katp+1)./(FBP./p.PYK1_Kf16bp+1)).^...
    p.PYK1_hill+(PEP./p.PYK1_Kpep+1).^p.PYK1_hill)));
v_LDH = p.LDH1_Vm .*(PYR .* NADH - LAC .* NAD ./ p.LDH1_Keq);

% mass balances (assumed ideally mixed batch system configuration)
vPYK(1) = - v_PYK;
vPYK(2) = - v_LDH;
vPYK(3) = 0;
vPYK(4) = - v_PYK;
vPYK(5) = + v_PYK;
vPYK(6) = + v_LDH;
vPYK(7) = + v_PYK - v_LDH;
vPYK(8) = + v_LDH;

vPYK = vPYK';
end