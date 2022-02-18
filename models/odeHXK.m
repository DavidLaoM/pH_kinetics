function vHXK = odeHXK(tspan,y0,p,f,pHdata,setup)

% select initial points
GLCi = y0(1); % mM
G6P = y0(2); % mM
PG6 = y0(3); % mM
ATP = y0(4); % mM
ADP = y0(5); % mM
NADP = y0(6); % mM
NADPH = y0(7); % mM

% calulate v (rateEquations)
v_GLK = (p.HXK1_Vm .* (ATP.*GLCi-((ADP.*G6P)./p.HXK1_Keq))) ./ ...
    ((p.HXK1_Katp.*p.HXK1_Kglc) .* (1+ATP./p.HXK1_Katp+ADP./p.HXK1_Kadp) .* ...
    (1+GLCi./p.HXK1_Kglc+G6P./p.HXK1_Kg6p+0./p.HXK1_Kt6p));
v_G6PDH = p.G6PDH_Vm .* (NADP .* G6P - NADPH .* PG6 ./ p.G6PDH_Keq);

% mass balances (assumed ideally mixed batch system configuration)
vHXK(1) = - v_GLK;
vHXK(2) = + v_GLK - v_G6PDH;
vHXK(3) = + v_G6PDH;
vHXK(4) = - v_GLK;
vHXK(5) = + v_GLK;
vHXK(6) = - v_G6PDH;
vHXK(7) = + v_G6PDH;

vHXK=vHXK';
end
