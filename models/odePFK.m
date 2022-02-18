function vPFK = odePFK(tspan,y0,p,f,pHdata,setup)
% disp(tspan);
% if tspan >= 360
%     disp('uno');
% end

% select initial points
FBP = y0(1); F16BP = FBP;% mM
DHAP = y0(2); % mM
GAP = y0(3); % mM
G3P = y0(4); % mM
NADH = y0(5); % mM
NAD = y0(6); % mM
F6P = y0(7); % mM
ATP = y0(8); % mM
ADP = y0(9); % mM
F26BP = 0.1; % mM
AMP = 0;

% calulate v (rateEquations)
for casePFK = 1
%     p.PFK_Kf6p = 0.12821;
%     p.PFK_Kf6p = 0.001;
%     p.PFK_Katp = 0.001;
%     p.PFK_Kf12bp = 0.001;
%     p.PFK_Kadp = 0.001;
    PFK_nom=(p.PFK_Vm.*p.PFK_gR.*(F6P./p.PFK_Kf6p).*(ATP./p.PFK_Katp).*(1+(F6P./p.PFK_Kf6p)+(ATP./p.PFK_Katp)+p.PFK_gR.*((F6P./p.PFK_Kf6p).*(ATP./p.PFK_Katp))));
    PFK_denom=(1+F6P./p.PFK_Kf6p+ATP./p.PFK_Katp+(p.PFK_gR.*(F6P./p.PFK_Kf6p).*(ATP./p.PFK_Katp))).^2+...
        p.PFK_L.*...
        ((1+p.PFK_Ciatp.*(ATP./p.PFK_Kiatp))./(1+ATP./p.PFK_Kiatp)).^2.*...
        ((1+p.PFK_Camp.*(AMP./p.PFK_Kamp))./(1+AMP./p.PFK_Kamp)).^2.*...
        ((1+((p.PFK_Cf26bp*F26BP)./(p.PFK_Kf26bp))+((p.PFK_Cf16bp.*F16BP)./(p.PFK_Kf16bp)))./(1+(F26BP./p.PFK_Kf26bp)+(F16BP./p.PFK_Kf16bp))).^2.*...
        (1+p.PFK_Catp.*(ATP./p.PFK_Katp)).^2;
    v_PFK = PFK_nom ./ PFK_denom;
end
% % % % % % % % % % for casePFK = 2
% % % % % % % % % %     v_PFK = p.PFK_Vm .* (F6P .* ATP - (FBP .* ADP)./p.PFK_Keq);
% % % % % % % % % % end
% % % % for casePFK = 2
% % % % % % % %     v_PFK = p.PFK_Vm .* F6P / 16;
% % % %     v_PFK = p.PFK_Vm .* F6P;
% % % % end
% % % % % % % % % % for caseSimplified_PFK = 1
% % % % % % % % % %     v_PFK = p.PFK_Vm;
% % % % % % % % % % end
% for casePFK = 4
%     p.PFK_Kadp = 100; % average from only 2 available values in BRENDA database. From other organisms, and not yeast.
%     v_PFK = (p.PFK_Vm .* (ATP.*F6P-((ADP.*FBP)./p.PFK_Keq))) ./ ...
%         ((p.PFK_Katp.*p.PFK_Kf6p) .* (1+ATP./p.PFK_Katp+ADP./p.PFK_Kadp) .* ...
%         (1+F6P./p.PFK_Kf6p+FBP./p.PFK_Kf16bp+0./p.PFK_Kf26bp));
% %     v_PFK = (p.PFK_Vm .* (ATP.*F6P-((0)./p.PFK_Keq))) ./ ...
% %         ((p.PFK_Katp.*p.PFK_Kf6p) .* (1+ATP./p.PFK_Katp+0./p.PFK_Kadp) .* ...
% %         (1+F6P./p.PFK_Kf6p+0./p.PFK_Kf16bp+0./p.PFK_Kf26bp));
% end
v_ALD = p.FBA1_Vm .* (FBP-(GAP.*DHAP)./p.FBA1_Keq);
v_GPD = p.GPD1_Vm .* (DHAP .* NADH - G3P .* NAD ./ p.GPD1_Keq);
v_TPI = p.TPI1_Vm .* (DHAP - GAP ./ p.TPI1_Keq);

% mass balances (assumed ideally mixed batch system configuration)
vPFK(1) = - v_ALD + v_PFK;
vPFK(2) = + v_ALD - v_GPD - v_TPI;
vPFK(3) = + v_ALD + v_TPI;
vPFK(4) = + v_GPD;
vPFK(5) = - v_GPD;
vPFK(6) = + v_GPD;
vPFK(7) = - v_PFK;
vPFK(8) = - v_PFK;
vPFK(9) = + v_PFK;

vPFK=vPFK';
end