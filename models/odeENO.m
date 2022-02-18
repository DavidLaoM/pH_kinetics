function vENO = odeENO(tspan,y0,p,f,pHdata,setup)

% select initial points
P2G = y0(1);
PEP = y0(2);

% calulate v (rateEquations)
v_ENO=((p.ENO1_vm./p.ENO1_K2pg).*(P2G-PEP./p.ENO1_Keq))./...
    (1+P2G./p.ENO1_K2pg+PEP./p.ENO1_Kpep);

% mass balances (assumed ideally mixed batch system configuration)
vENO(1) = - v_ENO;
vENO(2) = + v_ENO;

vENO=vENO';
end
